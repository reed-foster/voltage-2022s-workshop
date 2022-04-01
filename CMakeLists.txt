cmake_minimum_required(VERSION 3.13)

include(pico_sdk_import.cmake)

project(adc_fft_project C CXX ASM)

set(PICO_CXX_ENABLE_EXCEPTIONS 1)
set(PICO_CXX_ENABLE_RTTI 1)

pico_sdk_init()

include("PicoLed/PicoLed.cmake")

add_executable(specan specan.cpp)
add_library(kiss_fftr kiss_fftr.c)
add_library(kiss_fft kiss_fft.c)
add_compile_definitions(FIXED_POINT=8)

target_link_libraries(kiss_fftr kiss_fft)

pico_enable_stdio_usb(specan 1)
pico_enable_stdio_uart(specan 1)


target_link_libraries(specan
	pico_stdlib
	hardware_adc
	hardware_dma
	kiss_fftr
	PicoLed
	)

pico_add_extra_outputs(specan)
