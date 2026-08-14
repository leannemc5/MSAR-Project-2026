import numpy as np
import matplotlib.pyplot as plt

# constants
h = 6.626e-34
c = 3.0e8
kB = 1.381e-23

# planck function
def planck(wavelength, T):

    exponent = (h*c / (wavelength*kB*T))

    intensity = ((2*h*c**2) / (wavelength**5 * (np.exp(exponent)-1)))

    return intensity

# wavelength grid
wavelengths_nm = np.linspace(300, 1000, 5000)
wavelengths_m = wavelengths_nm * 1e-9

# blackbody
T = 6000
spectrum = planck(wavelengths_m, T)
spectrum /= np.max(spectrum)

# rayleigh scattering
airmass = 2

rayleigh = (550 / wavelengths_nm)**4

rayleigh /= np.max(rayleigh)

k_lambda = 0.4 * rayleigh

transmission = np.exp(-k_lambda * airmass)

extincted_spectrum = spectrum * transmission

extincted_spectrum /= np.max(extincted_spectrum)

# plot
plt.figure(figsize=(10,6))
plt.plot(wavelengths_nm, spectrum, label="Original")
plt.plot(wavelengths_nm, extincted_spectrum, label="Rayleigh Extinction")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Intensity")
plt.title("Wavelength Dependent Atmospheric Extinction")
plt.legend()
plt.grid(True)
plt.show()
