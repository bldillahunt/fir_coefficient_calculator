import math
import numpy as np
import sys
import matplotlib.pyplot as plt
from numpy import fft
from scipy import signal

fft_data_size = 2097152

coefficient_count = 0

while (coefficient_count <= 0) or (coefficient_count > 256):
    coefficient_count = int(input('Enter the number of taps: '))

sample_rate = 250e+6

# Create the ideal frequency response
# High pass starting at 1 MHz

ideal_response = []

filter_type = ' '

while (filter_type != 'l') and (filter_type != 'h'):
    print('Enter l for lowpass or h for highpass')
    filter_type = input('Enter filter type: ')

cutoff_frequency = 0 

while (cutoff_frequency <= 0):
    cutoff_frequency = int(input('Enter cutoff frequency: '))

cutoff_count = int((cutoff_frequency * fft_data_size)/sample_rate)

enable_sweep = ' '

while (enable_sweep != 'n') and (enable_sweep != 'y'):
    print('Enable or disable frequency sweep... ')
    enable_sweep = input('Enter "y" to enable, "n" to disable: ')

for i in range(cutoff_count):
    if (filter_type == 'h'):
        ideal_response.append(0)
    else:
        ideal_response.append(1)

for i in range(cutoff_count+1, fft_data_size):
    if (filter_type == 'l'):
        ideal_response.append(0)
    else:
        ideal_response.append(1)

# Take the inverse FFT to get the impulse response
impulse_response = np.fft.ifft(ideal_response)

sampling_rate = 1
n = len(impulse_response)
frequencies = np.fft.fftfreq(n, d=1/sampling_rate)

fft_shifted = np.fft.fftshift(impulse_response)
frequencies_shifted = np.fft.fftshift(frequencies)
magnitude_spectrum = fft_shifted    #np.abs(fft_shifted)

spectrum_real = magnitude_spectrum.real
spectrum_imag = magnitude_spectrum.imag

plt.figure(figsize=(10, 5))
plt.plot(frequencies_shifted, spectrum_real, color='blue')
plt.plot(frequencies_shifted, spectrum_imag, color='red')
plt.title('Magnitude Spectrum of the Signal')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.grid(True)
plt.show()

# Create window function
coefficients = []
sample_midpoint = fft_data_size/2

# Determine if the number of coefficients is even or odd
if (coefficient_count%2 == 0):
# Even number of coefficients
    negative_start = int(sample_midpoint) - int(coefficient_count/2) + 1
    positive_stop = int(sample_midpoint) + int(coefficient_count/2)

    for i in range(negative_start, int(sample_midpoint)+1):
        coefficients.append(magnitude_spectrum[i])
    
    for i in range(int(sample_midpoint)+1, positive_stop+1):
        coefficients.append(magnitude_spectrum[i])

else:
# Odd number of coefficients
# Start with the negative side and increment
    negative_start = int(sample_midpoint) - int(coefficient_count/2)
    positive_stop = int(sample_midpoint) + int(coefficient_count/2)
    midpoint_index = int(sample_midpoint)

    for i in range(negative_start, int(sample_midpoint)):
        coefficients.append(magnitude_spectrum[i])

    coefficients.append(magnitude_spectrum[midpoint_index])

    for i in range(int(sample_midpoint), positive_stop):
        coefficients.append(magnitude_spectrum[i])

print('length of coefficients = ', len(coefficients))

# Put coefficients into a file
with open('fir_coefficients.txt', 'w') as f:
    for i in range(0, len(coefficients)):
        print(coefficients[i], sep=',', file=f)

# Generate the waveform
output_bits = 28        # this gives 1 Hz resolution
rom_depth = 2**output_bits
accumulator_bits = output_bits*2

print('Enter time domain signal parameters...')

if (enable_sweep == 'y'):
    sweep_start_freq = input('Enter start frequency: ')
    sweep_stop_freq = input('Enter stop frequency: ')
    default_freq = 0
else:
    sweep_start_freq = 0
    sweep_stop_freq = 0
    default_freq = input('Enter default frequency: ')

frequency_sweep_time = 50e-6
frequency_sweep_sample_rate = sample_rate
delta_theta_low = (int(sweep_start_freq)*2**output_bits)/frequency_sweep_sample_rate
delta_theta_high = (int(sweep_stop_freq)*2**output_bits)/frequency_sweep_sample_rate
delta_theta_default = (int(default_freq)*2**output_bits)/frequency_sweep_sample_rate
clocks_per_sample_time = frequency_sweep_time*frequency_sweep_sample_rate
step_count_float = clocks_per_sample_time/(delta_theta_high-delta_theta_low)
step_count = int(step_count_float)
frequency_step = delta_theta_high-delta_theta_low
frequency_increment = float(frequency_step)/clocks_per_sample_time

print('delta_theta_low = ', delta_theta_low)
print('delta_theta_high = ', delta_theta_high)
print('delta_theta_default = ', delta_theta_default)
print('clocks_per_sample_time = ', clocks_per_sample_time)
print('step_count = ', step_count)
print('frequency_step = ', frequency_step)
print('frequency_increment = ', frequency_increment)

# First, create the lookup table
amplitude = 2**(output_bits - 1) - 1 # For signed output
#  t       S
# ---- = ----   => S = (2*PI*t)/2^16
# 2^16   2*PI
#sin_rom = [int(amplitude * np.sin(2 * np.pi * i / rom_depth)) for i in range(rom_depth)]
#cos_rom = [int(amplitude * np.cos(2 * np.pi * i / rom_depth)) for i in range(rom_depth)]
sin_rom = amplitude * np.sin(np.linspace(0, 2 * np.pi, rom_depth, endpoint=False))
cos_rom = amplitude * np.cos(np.linspace(0, 2 * np.pi, rom_depth, endpoint=False))

print('Calculated ROM tables')

# Second, create the phase accumulator
#ftw = int(delta_theta_default)  #int((frequency * (2**accumulator_bits)) / frequency_sweep_sample_rate)
phase_accumulator = 0

#ftw = int(delta_theta_default)
dds_output_i = []
dds_output_q = []
ftw_storage = []

# -------- Frequency sweep enabled --------
if (enable_sweep == 'y'):
    ftw = int(delta_theta_low*2**(accumulator_bits-output_bits))

    for i in range(0, int(clocks_per_sample_time)):
        phase_accumulator = (phase_accumulator + ftw) % float(2**accumulator_bits)
        rom_index = int(phase_accumulator) >> (accumulator_bits - int(np.log2(rom_depth)))
        dds_output_i.append(sin_rom[rom_index])
        dds_output_q.append(cos_rom[rom_index])
        ftw_storage.append(int(ftw)>>16)
        ftw = ftw + int(frequency_increment*2**(accumulator_bits-output_bits))
else:
# -------- Frequency sweep disabled --------
    ftw = int(delta_theta_default*2**(accumulator_bits-output_bits))

    for i in range(0, int(clocks_per_sample_time)):
        phase_accumulator = (phase_accumulator + ftw) % float(2**accumulator_bits)
        rom_index = int(phase_accumulator) >> (accumulator_bits - int(np.log2(rom_depth)))
        dds_output_i.append(sin_rom[rom_index])
        dds_output_q.append(cos_rom[rom_index])
        ftw_storage.append(int(ftw)>>16)

with open('dds_parameters.txt', 'w') as filename:
    for i in range(0, len(ftw_storage)):
        print(i, ftw_storage[i], file=filename)

print('DDS output length I = ', len(dds_output_i))
print('DDS output length Q = ', len(dds_output_q))

# Perform FFT of waveform before filter
complex_dds_data = np.array(dds_output_i) + 1j*np.array(dds_output_q)
pre_filter_waveform = np.fft.fft(complex_dds_data)

print(len(pre_filter_waveform))

# Obtain the frequencies of the filter
pre_filter_len = len(pre_filter_waveform)
pre_filter_rate = 1/frequency_sweep_sample_rate
pre_filter_freqs = np.fft.fftfreq(pre_filter_len, d=pre_filter_rate)

print('finished setting up pre-filter parameters')

# Shift the transform so that it is centered at 0
pre_filter_shifted = np.fft.fftshift(pre_filter_waveform)

# Shift the frequencies
pre_filter_freq_shifted = np.fft.fftshift(pre_filter_freqs)

# Apply the filter
a_coefficients = [1.0]
output_signal_i = signal.lfilter(coefficients, a_coefficients, dds_output_i)
output_signal_q = signal.lfilter(coefficients, a_coefficients, dds_output_q)

print('finished filtering data')

output_signal = output_signal_i + 1j*output_signal_q

dds_time = np.arange(0, len(complex_dds_data)/frequency_sweep_sample_rate, 1/frequency_sweep_sample_rate)

print('finished plot parameters for filtered time domain')

plt.plot(dds_time, complex_dds_data.real, color='blue')
plt.plot(dds_time, complex_dds_data.imag, color='red')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Unfiltered Time Domain Signal")
plt.grid(True)
plt.savefig("signal_time_domain.jpg")
plt.show()

plt.plot(dds_time, output_signal.real, color='blue')
plt.plot(dds_time, output_signal.imag, color='red')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Filtered Time Domain Signal")
plt.grid(True)
plt.savefig("filtered_time_domain.jpg")
plt.show()

print('finished time domain plot')

# Perform FFT of waveform after filter
post_filter_waveform = np.fft.fft(output_signal)

# Obtain the frequencies of the filter
post_filter_len = len(post_filter_waveform)
post_filter_rate = 1/frequency_sweep_sample_rate
post_filter_freqs = np.fft.fftfreq(post_filter_len, d=post_filter_rate)

# Shift the transform so that it is centered at 0
post_filter_shifted = np.fft.fftshift(post_filter_waveform)

# Shift the frequencies
post_filter_freq_shifted = np.fft.fftshift(post_filter_freqs)

pre_filter_db = 20 * np.log10(abs(pre_filter_shifted))

plt.figure(figsize=(10, 5))
plt.plot(pre_filter_freq_shifted, pre_filter_db, color='blue')
plt.title('FFT Output of Original Signal')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.grid(True)
plt.savefig("signal_fft.jpg")
plt.show()

post_filter_db = 20 * np.log10(abs(post_filter_shifted))

# Plot the input signal and output signal together
plt.figure(figsize=(10, 5))
#plt.plot(pre_filter_freq_shifted, pre_filter_shifted, color='blue')
plt.plot(post_filter_freq_shifted, post_filter_db, color='red')
plt.title('FFT Output of Filtered Signal')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.grid(True)
plt.savefig("filtered_signal_fft.jpg")
plt.show()

# Generate the frequency response of the filter
w, h = signal.freqz(coefficients, 1, whole=True)

magnitude_db = 20 * np.log10(abs(h))
phase_degrees = np.unwrap(np.angle(h)) * 180 / np.pi

plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.plot(w / (2 * np.pi), magnitude_db) # Plot against normalized frequency
plt.title('FIR Filter Frequency Response - Magnitude')
plt.xlabel('Normalized Frequency (cycles/sample)')
plt.ylabel('Magnitude (dB)')
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(w / (2 * np.pi), phase_degrees) # Plot against normalized frequency
plt.title('FIR Filter Frequency Response - Phase')
plt.xlabel('Normalized Frequency (cycles/sample)')
plt.ylabel('Phase (degrees)')
plt.grid()

plt.tight_layout()
plt.savefig("filter_response.jpg")
plt.show()

