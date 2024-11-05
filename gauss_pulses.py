import matplotlib.pyplot as plt
import numpy as np
import math

# Constants
c = 299792458  # Speed of light in m/s
h = 4.135667696e-15  # Planck's constant in eV·s
wavelength_center = 800e-9  # Central wavelength in meters
FWHM_intensity = 250e-15  # FWHM of intensity in seconds
line_thickness = 3  # Line thickness for plotting

period = 4000e-15
no_of_pulses = 1

sigma_t = FWHM_intensity / (2 * np.sqrt(2 * np.log(2)))  # Standard deviation of I(t)
t = np.linspace(-50000e-15, 50000e-15, 200000)  # Time array around pulse center

I_t = np.exp(-(t)**2 / (2 * sigma_t**2))
for pulse_no in range(1, int(math.floor(no_of_pulses/2)+1)):
    
    I_t += np.exp(-(t+period*pulse_no)**2 / (2 * sigma_t**2))
    I_t += np.exp(-(t-period*pulse_no)**2 / (2 * sigma_t**2))


f = np.fft.fftfreq(t.size, d=(t[1] - t[0]))  # Frequency array
I_f = np.fft.fft(I_t)  # Fourier transform of I(t)
f = np.fft.fftshift(f)  # Shift zero frequency to center
I_f = np.fft.fftshift(I_f)  # Shift corresponding I(f)

I_t /= np.max(np.abs(I_t))
I_f /= np.max(np.abs(2 * I_f))


plt.plot(t * 1e15, I_t, linewidth=line_thickness)  # Real part of E(t)
plt.figure()
plt.xlim([-5e12, 5e12])
plt.plot(f, np.abs(I_f), linewidth=line_thickness)  # Real part of E(t)
plt.show()

