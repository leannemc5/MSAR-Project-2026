import numpy as np
import matplotlib.pyplot as plt
from speclite import filters

# constants
h = 6.626e-34
c = 3.0e8
kB = 1.381e-23

# planck function
def planck(wavelength, T):

    exponent = (h * c / (wavelength * kB * T))

    intensity = (2 * h * c**2) / (wavelength**5 * (np.exp(exponent) - 1))

    return intensity

# wavelength grid
wavelengths_nm = np.linspace(300, 1000, 5000)
wavelengths_m = wavelengths_nm * 1e-9

# stellar spectrum
temperature = 6000

spectrum = planck(wavelengths_m, temperature)

# normalise 
spectrum = spectrum / np.max(spectrum)

# atmospheric extinction

airmass = 2.0

lambda0 = 550
k0 = 0.1

extinction_coeff = k0 * (lambda0 / wavelengths_nm)**4

atmospheric_transmission = np.exp(-extinction_coeff * airmass)

extincted_spectrum = spectrum * atmospheric_transmission

# filter response
bessell_v = filters.load_filter("bessell-V")

filter_response = np.interp(wavelengths_nm, bessell_v.wavelength / 10, bessell_v.response, left=0, right=0)

filtered_spectrum = extincted_spectrum * filter_response

# detector qe
qe_wavelength_nm = np.array([350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000])

qe_percent = np.array([20, 39, 55, 63, 72, 77, 85, 87, 86, 85, 80, 60, 38, 17])

qe_curve = np.interp(wavelengths_nm, qe_wavelength_nm, qe_percent / 100)

detected_spectrum = filtered_spectrum * qe_curve

# integrated flux
total_flux = np.trapezoid(detected_spectrum, wavelengths_m)

print(f"Integrated Flux = {total_flux:.6e}")

# simple photon count model
starting_photons = 100000

exposure_time = 10

expected_counts = total_flux * starting_photons * exposure_time

print(f"Expected Counts = {expected_counts:.0f}")

# atmospheric transparency
n_exposures = 500

cloud_factor = np.random.normal(loc=0.90, scale=0.03, size=n_exposures)

cloud_factor = np.clip(cloud_factor, 0, 1)

# light curve
light_curve = []

for transparency in cloud_factor:

    counts = expected_counts * transparency

    measured = np.random.poisson(counts)

    light_curve.append(measured)

light_curve = np.array(light_curve)

# transmission plot
plt.figure(figsize=(10,6))
plt.plot(wavelengths_nm, atmospheric_transmission, linewidth=2)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Transmission")
plt.title(f"Atmospheric Transmission (Airmass = {airmass})")
plt.grid(True)
plt.show()

# spectrum comparison
plt.figure(figsize=(10,6))
plt.plot(wavelengths_nm, spectrum, label="Original Spectrum")
plt.plot(wavelengths_nm, extincted_spectrum, label="After Atmospheric Extinction")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative Flux")
plt.title("Effect of Atmospheric Extinction")
plt.legend()
plt.grid(True)
plt.show()

# light curve
plt.figure(figsize=(10,6))
plt.plot(light_curve, marker=".", linewidth=0.8)

plt.xlabel("Exposure Number")
plt.ylabel("Detected Counts")
plt.title("Light Curve with Transparency Variations")
plt.grid(True)
plt.show()

# transparency
plt.figure(figsize=(10,6))
plt.plot(cloud_factor)

plt.xlabel("Exposure Number")
plt.ylabel("Transparency")
plt.title("Simulated Atmospheric Transparency")
plt.grid(True)
plt.show()

# statistics
print()
print(f"Mean Counts = {np.mean(light_curve):.2f}")
print(f"Standard Deviation = {np.std(light_curve):.2f}")
