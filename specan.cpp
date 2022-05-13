// Sample from the ADC continuously at a particular sample rate
// and then compute an FFT over the data, use power v. frequency to drive LEDs
//
// much of this code is from pico-examples/adc/dma_capture/dma_capture.c
// the rest is written by Alex Wulff (www.AlexWulff.com)

#include <stdio.h>
#include <math.h>
#include "pico/stdlib.h"
// multicore + semaphore
#include "pico/multicore.h"
#include "pico/sem.h"

// ADC/FFT
#include "hardware/adc.h"
#include "hardware/dma.h"
#include "kiss_fftr.h"
#include "_kiss_fft_guts.h"

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
// multiresolution analysis nfft scaling factors
#define HIGH_RES_FACTOR 8
#define MID_RES_FACTOR 4
// we'll offload the computation of the large FFT to a second core because it will take a while
// define SEMAPHOREDEBUG in order to print acquire/release messages from each core
//#define SEMAPHOREDEBUG
#endif

// FMAX/FMIN are maximum/minimum frequency displayed by LEDs
// make sure that MAX_BIN = FMAX/FSAMP*NSAMP and MIN_BIN = FMIN/FSAMP*NSAMP
#define MAX_BIN NSAMP/4
#define MIN_BIN NSAMP/480
#define FMAX 12000.0
#define FMIN 100.0

#define LED_STRIP_PIN 0

// globals
dma_channel_config cfg;
uint dma_chan;
// util for log-frequency spacing
double base;
double center_freqs[LED_STRIP_LENGTH];
uint16_t lower_bins[LED_STRIP_LENGTH];
uint16_t upper_bins[LED_STRIP_LENGTH];
// multiresolution analysis
#ifndef SHORT_WRGB
enum res_t { LR, MR, HR };
res_t resolution[LED_STRIP_LENGTH];
uint16_t hr_min_bin, hr_max_bin, hr_max_led_idx;
#endif /*SHORT_WRGB*/

// level curve data for fletcher-munson eq (TODO implement)
#ifdef SHORT_WRGB
uint16_t levelcurve[MAX_BIN];
#else /*SHORT_WRGB*/
uint16_t levelcurve_lr[MAX_BIN];
uint16_t levelcurve_mr[MID_RES_FACTOR*MAX_BIN];
uint16_t levelcurve_hr[HIGH_RES_FACTOR*MAX_BIN];
#endif /*SHORT_WRGB*/

#ifndef SHORT_WRGB
// shared global buffers and semaphore
kiss_fft_scalar * hr_fft_in;//[NSAMP*HIGH_RES_FACTOR];
kiss_fft_cpx * hr_fft_out;//[NSAMP*HIGH_RES_FACTOR];
semaphore_t fft_sem;
bool todo_fft;
#endif /*SHORT_WRGB*/

void setup();
void sample(uint8_t *capture_buf);
void compute_levelcurves();
double iso226(double f, double L_n);

#ifndef SHORT_WRGB
void core1_entry();

void core1_entry() {
    kiss_fftr_cfg hr_fft_cfg = kiss_fftr_alloc(NSAMP*HIGH_RES_FACTOR,false,0,0);
    // do large FFT on second core
    while (1) {
        bool do_fft = sem_acquire_timeout_us(&fft_sem, 100);
        if (do_fft) {
#ifdef SEMAPHOREDEBUG
            printf("1:acq\n");
#endif /*SEMAPHOREDEBUG*/
            if (todo_fft) {
                kiss_fftr(hr_fft_cfg, hr_fft_in, hr_fft_out);
                todo_fft = false;
            }
#ifdef SEMAPHOREDEBUG
            printf("1:rel\n");
#endif /*SEMAPHOREDEBUG*/
            sem_release(&fft_sem);
        }
        sleep_ms(2);
    }
    kiss_fft_free(hr_fft_cfg);
}
#endif /*SHORT_WRGB*/

int main() {
    uint8_t cap_buf[NSAMP];
    float led_psd[LED_STRIP_LENGTH];
  
    // track brightness of leds to give some persistence
    uint8_t led_brightness[LED_STRIP_LENGTH];
    // setup ports and outputs
    setup();
  
    // set up LED strip
#ifdef WRGB
    // color order for PicoLed::WRGB define isn't quite right; just do it manually
    auto ledStrip = PicoLed::addLeds<PicoLed::WS2812B>(pio0, 0, LED_STRIP_PIN, LED_STRIP_LENGTH, PicoLed::GREEN, PicoLed::RED, PicoLed::BLUE, PicoLed::WHITE);
#else /*WRGB*/
    auto ledStrip = PicoLed::addLeds<PicoLed::WS2812B>(pio0, 0, LED_STRIP_PIN, LED_STRIP_LENGTH, PicoLed::FORMAT_GRB);
#endif /*WRGB*/
    ledStrip.setBrightness(64);
    ledStrip.fill(PicoLed::RGB(255,0,0));
    ledStrip.show();

    sleep_ms(1000);
  
    uint8_t cycle = 0;
#ifndef SHORT_WRGB
    uint32_t last_sum[HIGH_RES_FACTOR] = {0};
#endif /*SHORT_WRGB*/

#ifdef SHORT_WRGB
    kiss_fft_scalar* fft_in = (kiss_fft_scalar *)malloc(NSAMP*sizeof(kiss_fft_scalar));
    kiss_fft_cpx* fft_out = (kiss_fft_cpx *)malloc(NSAMP*sizeof(kiss_fft_cpx));
    kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(NSAMP,false,0,0);
#else /*SHORT_WRGB*/
    kiss_fftr_cfg lr_fft_cfg = kiss_fftr_alloc(NSAMP,false,0,0);
    kiss_fftr_cfg mr_fft_cfg = kiss_fftr_alloc(NSAMP*MID_RES_FACTOR,false,0,0);
    // hr_fft_ibuf is used by core0
    // hr_fft_in is used by core1
    kiss_fft_scalar* hr_fft_ibuf = (kiss_fft_scalar *)malloc(NSAMP*HIGH_RES_FACTOR*sizeof(kiss_fft_scalar));
    hr_fft_in = (kiss_fft_scalar *)malloc(NSAMP*HIGH_RES_FACTOR*sizeof(kiss_fft_scalar));
    kiss_fft_scalar* mr_fft_in = (kiss_fft_scalar *)malloc(NSAMP*MID_RES_FACTOR*sizeof(kiss_fft_scalar));
    kiss_fft_scalar* lr_fft_in = (kiss_fft_scalar *)malloc(NSAMP*sizeof(kiss_fft_scalar));
    kiss_fft_cpx* mr_fft_out = (kiss_fft_cpx *)malloc(NSAMP*MID_RES_FACTOR*sizeof(kiss_fft_cpx));
    kiss_fft_cpx* lr_fft_out = (kiss_fft_cpx *)malloc(NSAMP*sizeof(kiss_fft_cpx));
    hr_fft_out = (kiss_fft_cpx *)malloc(NSAMP*HIGH_RES_FACTOR*sizeof(kiss_fft_cpx));
#endif /*SHORT_WRGB*/

    printf("MAX_BIN = %d, MIN_BIN = %d\n", MAX_BIN, MIN_BIN);
    printf("FMAX = %f, FMIN = %f\n", FMAX, FMIN);
    printf("levelcurves = [");
    for (int i = 0; i < MAX_BIN; i++) {
#ifdef SHORT_WRGB
        printf("%d", levelcurve[i]);
#else /*SHORT_WRGB*/
        printf("%d", levelcurve_lr[i]);
#endif /*SHORT_WRGB*/
        if (i < MAX_BIN - 1) {
            printf(", ");
            if (i % 32 == 31) { printf("\n"); }
        }
    }
    printf("]\n");
    
    printf("base = %f, FMAX/FMIN = %f\n", base, ((double)FMAX)/FMIN);
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
#ifdef SHORT_WRGB
            printf("(N/A,%d,%d) ", lower_bins[i*12+j], upper_bins[i*12+j]);
#else
            printf("(%d,%d,%d) ", resolution[i*12+j], lower_bins[i*12+j], upper_bins[i*12+j]);
#endif /*SHORT_WRGB*/
        }
        printf("\n");
    }
    printf("]\n");

#ifndef SHORT_WRGB
    // set up semaphore for sending data to and from core1 for computing the high res FFT
    todo_fft = true;
    sem_init(&fft_sem,1,1); // only one permit
    sleep_ms(100);
    // launch FFT engine on second core
    multicore_launch_core1(core1_entry);
    sleep_ms(100);
#endif /*SHORT_WRGB*/

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
#ifdef SHORT_WRGB
        for (int i = 0; i < NSAMP; i++) {
            fft_in[i] = cap_buf[i]-avg;
        }
        kiss_fftr(fft_cfg, fft_in, fft_out);
        // get LED brightness from PSD
        for (int i = 0; i < LED_STRIP_LENGTH; i++) {
            uint32_t led_power = 0;
            kiss_fft_cpx X;
            for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                X = fft_out[bin];
                led_power += levelcurve[bin]*(X.r*X.r + X.i*X.i);
            }
            led_psd[i] = 1e-3 * ((float)led_power / (upper_bins[i] - lower_bins[i] + 1));
        }
#else /*SHORT_WRGB*/
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
#ifdef SEMAPHOREDEBUG
            printf("0:acq\n");
#endif /*SEMAPHOREDEBUG*/
            if (!todo_fft) {
                // get data from FFT results and set LED brightness
                for (int i = 0; i < LED_STRIP_LENGTH; i++) {
                    if (resolution[i] != HR) { break; }
                    uint32_t led_power = 0;
                    for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                        kiss_fft_cpx X = hr_fft_out[bin];
                        led_power += levelcurve_hr[bin]*(X.r*X.r + X.i*X.i);
                    }
                    led_psd[i] = 1e-3 * (led_power / (upper_bins[i] - lower_bins[i] + 1));
                }

                // send core to do FFT
                memcpy(hr_fft_in, hr_fft_ibuf, sizeof(kiss_fft_scalar)*NSAMP*HIGH_RES_FACTOR);
                todo_fft = true;
            }
#ifdef SEMAPHOREDEBUG
            printf("0:rel\n");
#endif /*SEMAPHOREDEBUG*/
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
            kiss_fft_cpx X;
            if (resolution[i] == MR) {
                for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                    X = mr_fft_out[bin];
                    led_power += levelcurve_mr[bin]*(X.r*X.r + X.i*X.i);
                }
            } else {
                for (int bin = lower_bins[i]; bin <= upper_bins[i]; bin++) {
                    X = lr_fft_out[bin];
                    led_power += levelcurve_lr[bin]*(X.r*X.r + X.i*X.i);
                }
            }
            led_psd[i] = 1e-3 * ((float)led_power / (upper_bins[i] - lower_bins[i] + 1));
        }
#endif /*SHORT_WRGB*/
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
    }
#ifdef SHORT_WRGB
    kiss_fft_free(fft_cfg);  
    free(fft_in);  
    free(fft_out);  
#else /*SHORT_WRGB*/
    kiss_fft_free(lr_fft_cfg);
    kiss_fft_free(mr_fft_cfg);
    free(hr_fft_ibuf);
    free(hr_fft_in);
    free(mr_fft_in);
    free(lr_fft_in);
    free(mr_fft_out);
    free(lr_fft_out);
    free(hr_fft_out);
#endif /*SHORT_WRGB*/
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
    base = pow(((double)FMAX)/FMIN,1.0/((double)(LED_STRIP_LENGTH-1)));
#ifndef SHORT_WRGB
    hr_max_led_idx = 0;
    double lr_lim = (double)FSAMP*(1/((double)NSAMP));
    double mr_lim = (double)FSAMP*(1/((double)NSAMP*MID_RES_FACTOR));
#endif /*SHORT_WRGB*/
    for (int led_idx = 0; led_idx < LED_STRIP_LENGTH; led_idx++) {
        center_freqs[led_idx] = FMIN*pow(base, (double)led_idx);
#ifndef SHORT_WRGB
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
#endif /*SHORT_WRGB*/
    }

    for (int led_idx = 0; led_idx < LED_STRIP_LENGTH; led_idx++) {
        double bin_factor = ((double)NSAMP)/FSAMP;
#ifndef SHORT_WRGB
        if (resolution[led_idx] == HR) {
            bin_factor *= HIGH_RES_FACTOR;
        } else if (resolution[led_idx] == MR) {
            bin_factor *= MID_RES_FACTOR;
        }
#endif /*SHORT_WRGB*/
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
#ifndef SHORT_WRGB
    hr_min_bin = lower_bins[0];
    hr_max_bin = upper_bins[hr_max_led_idx];
#endif /*SHORT_WRGB*/

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
    // compensate for transfer function of speaker -> air -> microphone and weight based on perceived loudness
    // measured microphone response
    // f is frequency in Hz of a sinewave, a is amplitude of amplifier output in mVpp
    // f = [150 220 300 400 500 550 600 650 700 800 900 1000 1500 2000 2500 2800 3000 3500 4000 5000 7500 9000 12000 15000];
    // a = [100 200 330 600 750 710 690 640 640 640 700  660  570  650  460 1050 1050  640  250  660  190  190    90    90];
    // poles at 500Hz (bin 5) and 3500Hz (bin 37) for 512-point FFT
    // remeasured with good speaker and low freqency pole is gone
    // high frequency pole is probably still around 3500Hz
    //
    double L_n = 40; // assume 40 Phon for now, correct later with some sort of avg. power
    double L_p1k = iso226(1000, L_n);
    double f;
#ifdef SHORT_WRGB
    for (uint16_t i = 0 ; i < MAX_BIN; i++) {
        // pole at 5kHz
        f = i*FSAMP/NSAMP;
        levelcurve[i] = 256*(1+pow(f/3500,2))*pow(10, (L_p1k - iso226(f, L_n))/10);
    }
#else /*SHORT_WRGB*/
    for (uint16_t i = 0 ; i < HIGH_RES_FACTOR*MAX_BIN; i++) {
        if (i < MAX_BIN) {
            f = i*FSAMP/NSAMP;
            levelcurve_lr[i] = 256*(1+pow(f/3500,2))*pow(10, (L_p1k - iso226(f, L_n))/10);
        }
        if (i < MID_RES_FACTOR*MAX_BIN) {
            f = i*FSAMP/(NSAMP*MID_RES_FACTOR);
            levelcurve_mr[i] = 256*(1+pow(f/3500,2))*pow(10, (L_p1k - iso226(f, L_n))/10);
        }
        f = i*FSAMP/(NSAMP*HIGH_RES_FACTOR);
        levelcurve_hr[i] = 256*(1+pow(f/3500,2))*pow(10, (L_p1k - iso226(f, L_n))/10);
#endif /*SHORT_WRGB*/
    }
}

double iso226(double f, double L_n) {
    // iso226 loudness calculation
    double f0[29] = {20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500};
    double a_f0[29] = {0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349, 0.330, 0.315, 0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250, 0.246, 0.244, 0.243, 0.243, 0.243, 0.242, 0.242, 0.245, 0.254, 0.271, 0.301};
    double L_u0[29] = {-31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3, -8.1, -6.2, -4.5, -3.1, -2.0, -1.1, -0.4, 0.0, 0.3, 0.5, 0.0, -2.7, -4.1, -1.0, 1.7, 2.5, 1.2, -2.1, -7.1, -11.2, -10.7, -3.1};
    double T_f0[29] = {78.5, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9, 14.4, 11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3, -4.2, -6.0, -5.4, -1.5, 6.0, 12.6, 13.9, 12.3};
    double a_f, L_u, T_f;
    // linear interpolation of parameters
    int f0_low;
    for (f0_low = 0; f0[f0_low+1] <= f && f0_low < 28; f0_low++) { }
    a_f = a_f0[f0_low] + (a_f0[f0_low+1]-a_f0[f0_low])/(f0[f0_low+1]-f0[f0_low])*(f-f0[f0_low]);
    L_u = L_u0[f0_low] + (L_u0[f0_low+1]-L_u0[f0_low])/(f0[f0_low+1]-f0[f0_low])*(f-f0[f0_low]);
    T_f = T_f0[f0_low] + (T_f0[f0_low+1]-T_f0[f0_low])/(f0[f0_low+1]-f0[f0_low])*(f-f0[f0_low]);
    double A_F = 0.00447 * (pow(10, 0.025*L_n) - 1.15) + pow((0.4 * pow(10, ((T_f + L_u) / 10) - 9)), a_f);
    double L_p = (10 / a_f) * log10(A_F) - L_u + 94;
    return L_p;
}

