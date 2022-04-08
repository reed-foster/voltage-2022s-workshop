# Voltage 2022S CPW Workshop

Code adapted from [ForsakenNGS's PicoLED library](https://github.com/ForsakenNGS/PicoLED) and [AlexFWulff's pico-playground](https://github.com/AlexFWulff/awulff-pico-playground).

## What does this do?
This code turns your Pi Pico into a crude audio spectrum analyzer.
A microphone should be connect to one of the ADC channels (default channel 0) and a NeoPixel LED strip should be connected to a GPIO pin (default pin 0).
The Pi Pico uses DMA to stream the ADC data into memory, then computes a 512-point, 8-bit fixed point real [Fast Fourier Transform (FFT)](https://en.wikipedia.org/wiki/Fast_Fourier_transform) using the kiss FFT library.
Then, the 256 FFT bins (corresponding to frequencies from 0Hz up to Fsamp/2) are used to modulate the brightness and hue of 16 (or more, configurable through preprocessor define #LED_STRIP_LENGTH) LEDs.
The LEDs are mapped logarithmically to the range of frequencies available between 0Hz and Fsamp/2.
That is, a note which lights up an LED with index `i` played an octave higher would light up LED `i+1`, played two octaves higher would light up LED `i+2`, and so on.
The brightness and color of each LED is determined by the average [Power Spectral Density (PSD)](https://en.wikipedia.org/wiki/Spectral_density) over the frequency range assigned to that LED.
## Setup

![Breadboard](breadboard.png)

Ensure that you've installed the Pi Pico SDK. Follow the guide [here](https://datasheets.raspberrypi.com/pico/getting-started-with-pico.pdf)

Symlink the PicoLed library into the repo root
Add the lines to your `.bashrc`:
```
export PICO_SDK_PATH=/path/to/pico/pico-sdk
export PICO_EXAMPLES_PATH=/path/to/pico/pico-examples
```
replacing `/path/to/pico` with the path to the directory you installed the sdk in.

### Building a binary
```
$ mkdir build && cd !$
$ cmake ..
$ make
```

### Uploading the binary
```
# fdisk -l
...
Device     Boot Start    End Sectors  Size Id Type
/dev/sda1           1 262143  262143  128M  e W95 FAT16 (LBA)
# mount /dev/sda1 /mnt/pico
# pwd
~/voltage-2022s-workshop/build
# cp specan.uf2 /mnt/pico
```
