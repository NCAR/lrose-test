#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Sampling parameters
# ---------------------------------------------------------

fs = 1000.0
n_samples = 1024

t = np.arange(n_samples) / fs

# ---------------------------------------------------------
# Four complex sinusoids
# ---------------------------------------------------------

freqs = [50.0, 120.0, -200.0, 310.0]
amps  = [1.0, 0.7, 0.5, 0.3]

x = np.zeros(n_samples, dtype=complex)

for freq, amp in zip(freqs, amps):
    x += amp * np.exp(1j * 2.0 * np.pi * freq * t)

# ---------------------------------------------------------
# Complex FFT
# ---------------------------------------------------------

X = np.fft.fft(x)
f = np.fft.fftfreq(n_samples, d=1.0 / fs)

# Put negative frequencies on the left,
# zero in the middle, positive frequencies on the right.
X_shift = np.fft.fftshift(X)
f_shift = np.fft.fftshift(f)

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

fig, axes = plt.subplots(4, 1, figsize=(11, 10))

# Show only the first 0.1 seconds of the time series,
# otherwise the oscillations become visually crowded.
time_mask = t <= 0.1

axes[0].plot(t[time_mask], x.real[time_mask])
axes[0].set_ylabel("Real")
axes[0].set_title("Complex time series")
axes[0].grid(True)

axes[1].plot(t[time_mask], x.imag[time_mask])
axes[1].set_ylabel("Imaginary")
axes[1].set_xlabel("Time (s)")
axes[1].grid(True)

axes[2].plot(f_shift, X_shift.real)
axes[2].set_ylabel("Real")
axes[2].set_title("Complex FFT")
axes[2].grid(True)

axes[3].plot(f_shift, X_shift.imag)
axes[3].set_ylabel("Imaginary")
axes[3].set_xlabel("Frequency (Hz)")
axes[3].grid(True)

plt.tight_layout()
plt.show()

