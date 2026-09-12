#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Sampling parameters
# ---------------------------------------------------------

fs = 1000.0       # sample rate, Hz
n_samples = 1024

t = np.arange(n_samples) / fs

# ---------------------------------------------------------
# Four complex waveforms
# ---------------------------------------------------------

freqs = [50.0, 120.0, -200.0, 310.0]   # Hz
amps  = [1.0, 0.7, 0.5, 0.3]

x = np.zeros(n_samples, dtype=complex)

for freq, amp in zip(freqs, amps):
    x += amp * np.exp(1j * 2.0 * np.pi * freq * t)

# ---------------------------------------------------------
# Complex FFT
# ---------------------------------------------------------

X = np.fft.fft(x)

# Frequency corresponding to each FFT bin
f = np.fft.fftfreq(n_samples, d=1.0/fs)

# Shift negative frequencies to the left, positive to right
X = np.fft.fftshift(X)
f = np.fft.fftshift(f)

# Power spectrum
power = np.abs(X)**2

# Normalize for convenience
power /= power.max()

# Convert to dB
power_db = 10.0 * np.log10(np.maximum(power, 1.0e-12))

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

plt.figure(figsize=(10, 5))

plt.plot(f, power_db)

plt.xlabel("Frequency (Hz)")
plt.ylabel("Power (dB)")
plt.title("Spectrum of complex time series")

plt.grid(True)
plt.ylim(-80, 5)

plt.show()
