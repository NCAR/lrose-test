#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Radar parameters
# ---------------------------------------------------------

wavelength = 0.10
prt = 0.001
prf = 1.0 / prt

n_samples = 128
t = np.arange(n_samples) * prt

rng = np.random.default_rng(1)

# ---------------------------------------------------------
# Weather parameters
# ---------------------------------------------------------

weather_mean_velocity = 10.0    # m/s
weather_width = 2.0             # m/s, standard deviation
weather_power = 1.0

n_weather = 200

# Draw scatterer velocities
weather_velocity = rng.normal(
    weather_mean_velocity,
    weather_width,
    n_weather
)

weather_fd = 2.0 * weather_velocity / wavelength

# Random initial phases
weather_phase = rng.uniform(
    0.0, 2.0 * np.pi, n_weather
)

# Equal-power scatterers
weather_amp = np.sqrt(weather_power / n_weather)

weather = np.zeros(n_samples, dtype=complex)

for fd, phase in zip(weather_fd, weather_phase):
    weather += weather_amp * np.exp(
        1j * (2.0 * np.pi * fd * t + phase)
    )

# ---------------------------------------------------------
# Clutter parameters
# ---------------------------------------------------------

clutter_mean_velocity = 0.0
clutter_width = 0.25             # m/s
clutter_power = 100.0            # 20 dB above weather

n_clutter = 200

clutter_velocity = rng.normal(
    clutter_mean_velocity,
    clutter_width,
    n_clutter
)

clutter_fd = 2.0 * clutter_velocity / wavelength

clutter_phase = rng.uniform(
    0.0, 2.0 * np.pi, n_clutter
)

clutter_amp = np.sqrt(clutter_power / n_clutter)

clutter = np.zeros(n_samples, dtype=complex)

for fd, phase in zip(clutter_fd, clutter_phase):
    clutter += clutter_amp * np.exp(
        1j * (2.0 * np.pi * fd * t + phase)
    )

# ---------------------------------------------------------
# Combined signal
# ---------------------------------------------------------

x = weather + clutter

# ---------------------------------------------------------
# FFT
# ---------------------------------------------------------

X = np.fft.fft(x)
freq = np.fft.fftfreq(n_samples, d=prt)

X = np.fft.fftshift(X)
freq = np.fft.fftshift(freq)

velocity = freq * wavelength / 2.0

power = np.abs(X)**2
power_db = 10.0 * np.log10(
    np.maximum(power / power.max(), 1.0e-12)
)

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

