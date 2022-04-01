// Sample from the ADC continuously at a particular sample rate
// and then compute an FFT over the data, use power v. frequency to drive LEDs
//
// much of this code is from pico-examples/adc/dma_capture/dma_capture.c
// the rest is written by Alex Wulff (www.AlexWulff.com)

#include <stdio.h>
#include <math.h>
#include "pico/stdlib.h"

// ADC/multicore/FFT
#include "hardware/adc.h"
#include "hardware/dma.h"
#include "kiss_fftr.h"

// Neopixel LED Driver
#include <PicoLed.hpp>

// set this to determine sample rate
// 0     = 500,000 Hz
// 960   = 50,000 Hz
// 9600  = 5,000 Hz
#define CLOCK_DIV 1000
#define FSAMP 48000

// Channel 0 is GPIO26
#define CAPTURE_CHANNEL 0
#define LED_PIN 25

// BE CAREFUL: anything over about 9000 here will cause things
// to silently break. The code will compile and upload, but due
// to memory issues nothing will work properly
#define NSAMP 512

#define LED_STRIP_LENGTH 10
#define LED_STRIP_PIN 0

// globals
dma_channel_config cfg;
uint dma_chan;
// get fft bin index ranges for each LED
uint16_t low_bins[LED_STRIP_LENGTH];
uint16_t high_bins[LED_STRIP_LENGTH];
int8_t levelcurve[NSAMP/2];

void setup();
void sample(uint8_t *capture_buf);
float get_bin_log_base(uint16_t led_count, uint16_t fft_length);
void compute_levelcurve(int8_t *levelcurve);

int main() {
  uint8_t cap_buf[NSAMP];
  uint32_t power[NSAMP/2];
  uint64_t total_power = 0;
  uint64_t max_power = 0;
  kiss_fft_scalar fft_in[NSAMP]; // kiss_fft_scalar is a float
  kiss_fft_cpx fft_out[NSAMP];
  kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(NSAMP,false,0,0);

  // track brightness of leds to give some persistence
  uint8_t led_brightness[LED_STRIP_LENGTH];
  // setup ports and outputs
  setup();

  // set up LED strip
  auto ledStrip = PicoLed::addLeds<PicoLed::WS2812B>(pio0, 0, LED_STRIP_PIN, LED_STRIP_LENGTH, PicoLed::FORMAT_GRB);
  ledStrip.setBrightness(64);
  ledStrip.fill(PicoLed::RGB(255,0,0));
  ledStrip.show();

  while (1) {
    // get NSAMP samples at FSAMP
    sample(cap_buf);
    // fill fourier transform input while subtracting DC component
    uint32_t sum = 0;
    for (int i = 0; i < NSAMP; i++) {
        uint16_t val = cap_buf[i];
        sum += val;
    }
    // don't worry about rounding error; there's a bunch of noise already
    uint8_t avg = (uint8_t)(sum/NSAMP);
    for (int i = 0; i < NSAMP; i++) {
        fft_in[i] = cap_buf[i]-avg;
    }

    // compute fast fourier transform and get power for each bin
    // since fft_in[n] is real, we only need half the spectrum since it will be symmetric
    kiss_fftr(fft_cfg , fft_in, fft_out);
    for (int i = 0; i < NSAMP/2; i++) {
        power[i] = fft_out[i].r*fft_out[i].r+fft_out[i].i*fft_out[i].i;
        total_power += power[i];
    }
    if (total_power > max_power) {
        max_power = total_power;
    } else {
        total_power = total_power * 0.95;
    }
    // get LED brightness
    for (int i = 0; i < LED_STRIP_LENGTH; i++) {
        uint32_t power_led = 0;
        for (int bin = low_bins[i]; bin <= high_bins[i]; bin++) {
            power_led += power[bin]*levelcurve[bin];
        }
        // printf("led %02i power = %d\n", i, power_led);
        // convert to brightness
        uint8_t new_brightness;
        if (power_led >= 1500) {
            new_brightness = 255;
        } else {
            new_brightness = (uint8_t)((power_led*256)/1500);
        }
        if (new_brightness > led_brightness[i]) {
            led_brightness[i] = new_brightness;
        } else {
            led_brightness[i] = 0.9 * led_brightness[i];
        }
        new_brightness = led_brightness[i];
        // set brightness of led
        ledStrip.setPixelColor(i, PicoLed::RGB(new_brightness,new_brightness,new_brightness));
    }
    ledStrip.show();
  }

  // should never get here
  kiss_fft_free(fft_cfg);
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
  float base = get_bin_log_base(LED_STRIP_LENGTH, NSAMP/2);
  uint16_t count = 0;
  if (base) {
      for (int led_idx = 0; led_idx < LED_STRIP_LENGTH; led_idx++) {
          low_bins[led_idx] = count;
          printf("bin %d range = %5.2f - ", led_idx, (float)count*FSAMP/NSAMP);
          count += int(pow(base, led_idx) + 0.5);
          printf("%5.2f\n", (float)count*FSAMP/NSAMP);
          high_bins[led_idx] = count - 1;
      }
  } else {
      printf("error getting log base\n");
  }

  // precompute levelcurve data
  compute_levelcurve(levelcurve);

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

float get_bin_log_base(uint16_t led_count, uint16_t fft_length) {
    // determine what logarithmic spacing should be used
    // to fit fft_length frequency bins into led_count leds
    float increment = 0.1;
    // start with guess of base = 1
    // then increment base until sum(base^i, i = {0,led_count}) is as close to fft_length as possible
    for (float base = 1; base < fft_length; base += increment) {
        uint16_t count = 0;
        for (uint16_t led_idx = 0; led_idx < led_count; led_idx++) {
            count += int(pow(base, led_idx) + 0.5);
        }
        if (count > fft_length) {
            base -= increment;
            increment /= 10.0;
        } else if (count == fft_length) {
            return base;
        }
        if (increment < 0.000001) {
            return base - increment;
        }
    }
    return 0;
}

void compute_levelcurve(int8_t *levelcurve) {
    // measured microphone response
    // f is frequency in Hz of a sinewave, a is amplitude of amplifier output in mVpp
    // f = [150 220 300 400 500 550 600 650 700 800 900 1000 1500 2000 2500 2800 3000 3500 4000 5000 7500 9000 12000 15000];
    // a = [100 200 330 600 750 710 690 640 640 640 700  660  570  650  460 1050 1050  640  250  660  190  190    90    90];
    for (uint16_t i = 0 ; i < NSAMP/2; i++) {
        levelcurve[i] = 1;
    }
}

