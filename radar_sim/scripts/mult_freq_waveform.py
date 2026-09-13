#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Sampling parameters
# ---------------------------------------------------------

fs = 1000.0       # sample rate, Hz
n_samples = 1024
#n_samples = 112

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

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# ---------------------------------------------------------
# Plot time-series magnitude
# ---------------------------------------------------------

axes[0].plot(t, np.abs(x))
axes[0].set_xlabel("Time (s)")
axes[0].set_ylabel("Magnitude")
axes[0].set_title("Magnitude of complex time series")
axes[0].grid(True)

# ---------------------------------------------------------
# Plot spectrum
# ---------------------------------------------------------

axes[1].plot(f, power_db)
axes[1].set_xlabel("Frequency (Hz)")
axes[1].set_ylabel("Power (dB)")
axes[1].set_title("Spectrum of complex time series")
#axes[1].set_ylim(-180, 5)
axes[1].grid(True)

plt.tight_layout()
plt.show()

