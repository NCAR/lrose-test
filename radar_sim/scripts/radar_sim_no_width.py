#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Radar parameters
# ---------------------------------------------------------

wavelength = 0.10       # m
prt = 0.001             # s
prf = 1.0 / prt         # Hz

n_samples = 128

t = np.arange(n_samples) * prt

# ---------------------------------------------------------
# Weather
# ---------------------------------------------------------

weather_velocity = 10.0       # m/s, away from radar
weather_amp = 1.0

weather_fd = 2.0 * weather_velocity / wavelength

weather = weather_amp * np.exp(
    1j * 2.0 * np.pi * weather_fd * t
)

# ---------------------------------------------------------
# Ground clutter
# ---------------------------------------------------------

clutter_velocity = 0.0

# Clutter is 20 dB stronger in POWER.
# Therefore it is 10 times stronger in voltage amplitude.

clutter_db = 20.0
clutter_amp = weather_amp * 10.0**(clutter_db / 20.0)

clutter_fd = 2.0 * clutter_velocity / wavelength

clutter = clutter_amp * np.exp(
    1j * 2.0 * np.pi * clutter_fd * t
)

# ---------------------------------------------------------
# Combined received complex voltage
# ---------------------------------------------------------

x = weather + clutter

# ---------------------------------------------------------
# FFT
# ---------------------------------------------------------

X = np.fft.fft(x)

freq = np.fft.fftfreq(n_samples, d=prt)

X = np.fft.fftshift(X)
freq = np.fft.fftshift(freq)

power = np.abs(X)**2

# Normalize spectrum to maximum
power_db = 10.0 * np.log10(
    np.maximum(power / power.max(), 1.0e-12)
)

# Convert spectral frequency to radial velocity
velocity = freq * wavelength / 2.0

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# ---------------------------------------------------------
# Plot complex time-series magnitude
# ---------------------------------------------------------

axes[0].plot(t * 1000.0, np.abs(x))
axes[0].set_xlabel("Time (s)")
axes[0].set_ylabel("Voltage magnitude")
axes[0].set_title("Received complex time series")
axes[0].grid(True)

# ---------------------------------------------------------
# Plot spectrum
# ---------------------------------------------------------

axes[1].plot(freq, power_db)
axes[1].set_xlabel("Doppler frequency (Hz)")
axes[1].set_ylabel("Power (dB)")
axes[1].set_title("Doppler spectrum")
axes[1].grid(True)

plt.show()

print("PRF =", prf, "Hz")
print("Weather Doppler frequency =", weather_fd, "Hz")
print("Clutter Doppler frequency =", clutter_fd, "Hz")

print("Weather amplitude =", weather_amp)
print("Clutter amplitude =", clutter_amp)

