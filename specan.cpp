// Sample from the ADC continuously at a particular sample rate
// and then compute an FFT over the data, use power v. frequency to drive LEDs
//
// much of this code is from pico-examples/adc/dma_capture/dma_capture.c
// the rest is written by Alex Wulff (www.AlexWulff.com)

#include <stdio.h>
#include <math.h>
#include "pico/stdlib.h"
// multicore + queue
#include "pico/multicore.h"
#include "pico/util/queue.h"
#include "pico/sem.h"

// ADC/FFT
#include "hardware/adc.h"
#include "hardware/dma.h"
#include "kiss_fftr.h"

// Neopixel LED Driver
#include <PicoLed.hpp>

// colormaps
#include "colormaps.h"

// set this to determine sample rate
// 0     = 500,000 Hz
// 960   = 50,000 Hz
// 9600  = 5,000 Hz
#define CLOCK_DIV 1000
#define FSAMP 48000

// Channel 0 is GPIO26
#define CAPTURE_CHANNEL 0
#define LED_PIN 25

//#define SHORT_WRGB

// BE CAREFUL: anything over about 9000 for NSAMP will cause things
// to silently break. The code will compile and upload, but due
// to memory issues nothing will work properly
#ifdef SHORT_WRGB
#define NSAMP 512
#define LED_STRIP_LENGTH 16
#define WRGB 0
#else
#define NSAMP 512
#define LED_STRIP_LENGTH 120
// multiresolution analysis:
// HIGH_RES_FACTOR: X*NSAMP-point FFT for high resolution transform
// for NSAMP = 512, we'll have FFT data with a resolution of 96Hz/div
// MID_RES FFT will incrase that to 24Hz/div, and HIGH_RES FFT will have
// a resolution of 12Hz/div
#define HIGH_RES_FACTOR 16
#define MID_RES_FACTOR 4
// we'll offload the computation of the large FFT to a second core because it will take a while
#endif

// FMAX/FMIN are maximum/minimum frequency displayed by LEDs
#define FMAX 12000.0
#define FMIN 50
#define MAX_BIN NSAMP*FMAX/FSAMP
#define MIN_BIN NSAMP*FMIN/FSAMP

#define LED_STRIP_PIN 0

// globals
dma_channel_config cfg;
uint dma_chan;
// get fft bin index ranges for each LED
double base;
double center_freqs[LED_STRIP_LENGTH];
enum res_t { LR, MR, HR };
res_t resolution[LED_STRIP_LENGTH];
uint16_t lower_bins[LED_STRIP_LENGTH];
uint16_t upper_bins[LED_STRIP_LENGTH];
uint16_t hr_min_bin, hr_max_bin, hr_max_led_idx;
// level curve data for fletcher-munson eq (TODO implement)
uint16_t levelcurve_lr[NSAMP/2];
uint16_t levelcurve_mr[MID_RES_FACTOR*NSAMP/2];
uint16_t levelcurve_hr[HIGH_RES_FACTOR*NSAMP/2];

// shared global buffers and semaphore
kiss_fft_scalar hr_fft_in[NSAMP*HIGH_RES_FACTOR]; // kiss_fft_scalar is a float
kiss_fft_cpx hr_fft_out[NSAMP*HIGH_RES_FACTOR];
semaphore_t fft_sem;
bool todo_fft;

void setup();
void sample(uint8_t *capture_buf);
void compute_levelcurves();
void core1_entry();

void core1_entry() {
    kiss_fftr_cfg hr_fft_cfg = kiss_fftr_alloc(NSAMP*HIGH_RES_FACTOR,false,0,0);
    // do large FFT on second core
    while (1) {
        bool do_fft = sem_acquire_timeout_us(&fft_sem, 100);
        if (do_fft) {
            printf("1:acq\n");
            if (todo_fft) {
                kiss_fftr(hr_fft_cfg, hr_fft_in, hr_fft_out);
                todo_fft = false;
            }
            printf("1:rel\n");
            sem_release(&fft_sem);
        }
        sleep_ms(2);
    }
    kiss_fft_free(hr_fft_cfg);
}

int main() {
    uint8_t cap_buf[NSAMP];
    float led_psd[LED_STRIP_LENGTH];
    kiss_fft_scalar lr_fft_in[NSAMP]; // kiss_fft_scalar is a float
    kiss_fft_cpx lr_fft_out[NSAMP];
    kiss_fftr_cfg lr_fft_cfg = kiss_fftr_alloc(NSAMP,false,0,0);
    kiss_fft_scalar mr_fft_in[NSAMP*MID_RES_FACTOR]; // kiss_fft_scalar is a float
    kiss_fft_cpx mr_fft_out[NSAMP*MID_RES_FACTOR];
    kiss_fftr_cfg mr_fft_cfg = kiss_fftr_alloc(NSAMP*MID_RES_FACTOR,false,0,0);
    // keep a copy of the global input buffer for the highres FFT
    kiss_fft_scalar hr_fft_ibuf[NSAMP*HIGH_RES_FACTOR];
  
    // track brightness of leds to give some persistence
    uint8_t led_brightness[LED_STRIP_LENGTH];
    // setup ports and outputs
    setup();
  
    // set up LED strip
#ifdef WRGB
    // color order for PicoLed::WRGB define isn't quite right; just do it manually
    auto ledStrip = PicoLed::addLeds<PicoLed::WS2812B>(pio0, 0, LED_STRIP_PIN, LED_STRIP_LENGTH, PicoLed::GREEN, PicoLed::RED, PicoLed::BLUE, PicoLed::WHITE);
#else
    auto ledStrip = PicoLed::addLeds<PicoLed::WS2812B>(pio0, 0, LED_STRIP_PIN, LED_STRIP_LENGTH, PicoLed::FORMAT_GRB);
#endif
    ledStrip.setBrightness(64);
    ledStrip.fill(PicoLed::RGB(255,0,0));
    ledStrip.show();
  
    uint8_t cycle = 0;
    uint32_t last_sum[HIGH_RES_FACTOR] = {0};

    sleep_ms(1000);
    printf("base = %2.5f\n", base);
    printf("center_freqs = [\n");
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 12; j++) {
            printf("%f ", center_freqs[i*12+j]);
        }
        printf("\n");
    }
    printf("]\n");
    printf("(resolution, lower_bins, upper_bins) = [\n");
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 12; j++) {
            printf("(%d,%d,%d) ", resolution[i*12+j], lower_bins[i*12+j], upper_bins[i*12+j]);
        }
        printf("\n");
    }
    printf("]\n");
    
    // set up semaphore for sending data to and from core1 for computing the high res FFT
    todo_fft = true;
    sem_init(&fft_sem,1,1); // only one permit
    sleep_ms(100);
    // launch FFT engine on second core
    multicore_launch_core1(core1_entry);
    sleep_ms(100);

    while (1) {
        // get NSAMP samples at FSAMP
        sample(cap_buf);
        // fill fourier transform input while subtracting DC component
        uint32_t sum = 0;
        for (int i = 0; i < NSAMP; i++) {
            uint16_t val = cap_buf[i];
            sum += val;
        }
        // don't worry about non-symmetric rounding error; there's a bunch of noise already
        uint8_t avg = (uint8_t)(sum/NSAMP);
        for (int j = 0; j < HIGH_RES_FACTOR-1; j++) {
            last_sum[j+1] = last_sum[j];
        }
        last_sum[0] = sum;
        uint32_t mr_sum = 0;
        uint32_t hr_sum = 0;
        for (int j = 0; j < HIGH_RES_FACTOR; j++) {
            if (j < MID_RES_FACTOR) {
                mr_sum += last_sum[j];
            }
            hr_sum += last_sum[j];
        }
        uint8_t mr_avg = (uint8_t)(mr_sum/(NSAMP*MID_RES_FACTOR));
        uint8_t hr_avg = (uint8_t)(hr_sum/(NSAMP*HIGH_RES_FACTOR));
        for (int i = 0; i < NSAMP; i++) {
            lr_fft_in[i] = cap_buf[i]-avg;
            mr_fft_in[(cycle % MID_RES_FACTOR)*NSAMP+i] = cap_buf[i]-mr_avg;
            hr_fft_ibuf[(cycle % HIGH_RES_FACTOR)*NSAMP+i] = cap_buf[i]-hr_avg;
        }
        bool process_fft = sem_acquire_timeout_us(&fft_sem, 100);
        if (process_fft) {
            printf("0:acq\n");
            if (!todo_fft) {
                // get data from FFT results and set LED brightness
                printf("HR LED_PSD = [\n");
                for (int i = 0; i < LED_STRIP_LENGTH; i++) {
                    if (resolution[i] != HR) { break; }
                    uint32_t led_power = 0;
                    for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                        kiss_fft_cpx X = hr_fft_out[bin];
                        led_power += levelcurve_hr[bin]*(X.r*X.r + X.i*X.i);
                    }
                    led_psd[i] = 1e-3 * (led_power / (upper_bins[i] - lower_bins[i] + 1));
                    printf("%f, ", led_psd[i]);
                    if (i % 12 == 11) {
                        printf("\n");
                    }
                }
                printf("]\n");

                // send core to do FFT
                memcpy(hr_fft_in, hr_fft_ibuf, sizeof(kiss_fft_scalar)*NSAMP*HIGH_RES_FACTOR);
                todo_fft = true;
            }
            printf("0:rel\n");
            sem_release(&fft_sem);
        } else {
            sleep_ms(2);
        }
  
        // compute MID/LOW RES FFTs
        // since fft_in[n] is real, we only need half the spectrum since it will be symmetric
        kiss_fftr(lr_fft_cfg, lr_fft_in, lr_fft_out);
        kiss_fftr(mr_fft_cfg, mr_fft_in, mr_fft_out);

        // get LED brightness from PSD for non-highres PSD
        for (int i = hr_max_led_idx+1; i < LED_STRIP_LENGTH; i++) {
            uint32_t led_power = 0;
            uint8_t factor; // how much to multiply PSD by to account for frequency spacing of bins
            kiss_fft_cpx X;
            if (resolution[i] == MR) {
                factor = MID_RES_FACTOR;
                for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                    X = mr_fft_out[bin];
                    led_power += levelcurve_mr[bin]*(X.r*X.r + X.i*X.i);
                }
            } else {
                factor = 1;
                for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                    X = lr_fft_out[bin];
                    led_power += levelcurve_lr[bin]*(X.r*X.r + X.i*X.i);
                }
            }
            led_psd[i] = 1e-3 * ((float)led_power / (upper_bins[i] - lower_bins[i] + 1));
        }
        for (int i = 0; i < LED_STRIP_LENGTH; i++) {
            // convert to brightness
            float log_psd = log2(led_psd[i]);
            uint8_t color_idx;
            float new_brightness = pow(led_psd[i],2);
            if (new_brightness >= 10) {
                new_brightness = 10;
            }
            if (log_psd > 0) {
                color_idx = 255;
            } else if (log_psd > -8) {
                color_idx = ((log_psd+8)*255)/8;
            } else {
                color_idx = 0;
            }
            new_brightness = (255*(new_brightness+0.1))/10;
            if (new_brightness > led_brightness[i]) {
                led_brightness[i] = (uint8_t)new_brightness;
            } else {
                led_brightness[i] = 0.9 * led_brightness[i];
            }
            uint8_t b = led_brightness[i];
#ifdef WRGB
            b = b >> 4;
            b = b <= 1 ? 1 : b;
            uint8_t cb = (b <= 4 ? 4 : b) << 4;
            PicoLed::Color color = PicoLed::RGBW((cb*turbo[3*color_idx])/255, (cb*turbo[3*color_idx+1])/255, (cb*turbo[3*color_idx+2])/255, b);
#else
            b = b <= 32 ? 32 : b;
            PicoLed::Color color = PicoLed::RGB((b*turbo[3*color_idx])/255, (b*turbo[3*color_idx+1])/255, (b*turbo[3*color_idx+2])/255);
#endif
            // set brightness of led
            ledStrip.setPixelColor(i, color);
        }
        ledStrip.show();
        cycle++;
        cycle = cycle % MID_RES_FACTOR;
    }
  
    // should never get here
    kiss_fft_free(lr_fft_cfg);
    kiss_fft_free(mr_fft_cfg);
}

void sample(uint8_t *capture_buf) {
    adc_fifo_drain();
    adc_run(false);
        
    dma_channel_configure(dma_chan, &cfg,
    		capture_buf,    // dst
    		&adc_hw->fifo,  // src
    		NSAMP,          // transfer count
    		true            // start immediately
    		);

    gpio_put(LED_PIN, 1);
    adc_run(true);
    dma_channel_wait_for_finish_blocking(dma_chan);
    gpio_put(LED_PIN, 0);
}

void setup() {
    stdio_init_all();

    // calculate allocation of frequency ranges to leds
    base = pow((double)FMAX/FMIN,1.0/((double)(LED_STRIP_LENGTH-1)));
    hr_max_led_idx = 0;
    double lr_lim = (double)FSAMP*(0.5/((double)NSAMP));
    double mr_lim = (double)FSAMP*(0.5/((double)NSAMP*MID_RES_FACTOR));
    for (int led_idx = 0; led_idx < LED_STRIP_LENGTH; led_idx++) {
        center_freqs[led_idx] = FMIN*pow(base, (double)led_idx);
        if (led_idx == 0) {
            resolution[led_idx] = HR;
        } else {
            if ((center_freqs[led_idx] - center_freqs[led_idx-1]) > lr_lim) {
                resolution[led_idx] = LR;
            } else if ((center_freqs[led_idx] - center_freqs[led_idx-1]) > mr_lim) {
                resolution[led_idx] = MR;
            } else {
                resolution[led_idx] = HR;
                hr_max_led_idx = led_idx;
            }
        }
    }
    for (int led_idx = 0; led_idx < LED_STRIP_LENGTH; led_idx++) {
        double bin_factor = ((double)NSAMP)/FSAMP;
        if (resolution[led_idx] == HR) {
            bin_factor *= HIGH_RES_FACTOR;
        } else if (resolution[led_idx] == MR) {
            bin_factor *= MID_RES_FACTOR;
        }
        if (led_idx == 0) {
            lower_bins[led_idx] = (uint16_t)round(center_freqs[led_idx]*bin_factor);
        } else {
            lower_bins[led_idx] = (uint16_t)round(0.5*(center_freqs[led_idx] + center_freqs[led_idx-1])*bin_factor);
        }
        if (led_idx == LED_STRIP_LENGTH - 1) {
            upper_bins[led_idx] = (uint16_t)round(center_freqs[led_idx]*bin_factor);
        } else {
            upper_bins[led_idx] = (uint16_t)round(0.5*(center_freqs[led_idx] + center_freqs[led_idx+1])*bin_factor);
        }
    }
    hr_min_bin = lower_bins[0];
    hr_max_bin = upper_bins[hr_max_led_idx];

    // precompute levelcurve data
    compute_levelcurves();

    gpio_init(LED_PIN);
    gpio_set_dir(LED_PIN, GPIO_OUT);

    adc_gpio_init(26 + CAPTURE_CHANNEL);

    adc_init();
    adc_select_input(CAPTURE_CHANNEL);
    adc_fifo_setup(
	  	 true,    // Write each completed conversion to the sample FIFO
	  	 true,    // Enable DMA data request (DREQ)
	  	 1,       // DREQ (and IRQ) asserted when at least 1 sample present
	  	 false,   // We won't see the ERR bit because of 8 bit reads; disable.
	  	 true     // Shift each sample to 8 bits when pushing to FIFO
	  	 );

    // set sample rate
    adc_set_clkdiv(CLOCK_DIV);

    sleep_ms(1000);
    // Set up the DMA to start transferring data as soon as it appears in FIFO
    uint dma_chan = dma_claim_unused_channel(true);
    cfg = dma_channel_get_default_config(dma_chan);

    // Reading from constant address, writing to incrementing byte addresses
    channel_config_set_transfer_data_size(&cfg, DMA_SIZE_8);
    channel_config_set_read_increment(&cfg, false);
    channel_config_set_write_increment(&cfg, true);

    // Pace transfers based on availability of ADC samples
    channel_config_set_dreq(&cfg, DREQ_ADC);
}

void compute_levelcurves() {
    // measured microphone response
    // f is frequency in Hz of a sinewave, a is amplitude of amplifier output in mVpp
    // f = [150 220 300 400 500 550 600 650 700 800 900 1000 1500 2000 2500 2800 3000 3500 4000 5000 7500 9000 12000 15000];
    // a = [100 200 330 600 750 710 690 640 640 640 700  660  570  650  460 1050 1050  640  250  660  190  190    90    90];
    // poles at 500Hz (bin 5) and 3500Hz (bin 37) for 512-point FFT
    // remeasured with good speaker and low freqency pole is gone
    // high frequency pole is probably still around 3500Hz
    for (uint16_t i = 0 ; i < NSAMP/2; i++) {
        // pole at 5kHz
        levelcurve_lr[i] = 256*(1+pow((i*FSAMP/NSAMP)/3500,2));
        levelcurve_lr[i] = 256;
        levelcurve_mr[i] = 256;
        levelcurve_hr[i] = 256;
    }
}

