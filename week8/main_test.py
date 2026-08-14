import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import h, c, k
from scipy.interpolate import interp1d
from speclite import filters

# planck function
def planck(wavelength, T):

    exponent = (h * c) / (wavelength * k * T)

    intensity = ((2 * h * c**2) / (wavelength**5 * (np.exp(exponent) - 1)))

    return intensity

# wavelength grid
wavelengths_nm = np.linspace(300, 1100, 5000)
wavelengths_m = wavelengths_nm * 1e-9

# atmospheric model
airmass = 0.5

rayleigh = (550 / wavelengths_nm)**4

rayleigh /= np.max(rayleigh)

k_lambda = 0.4 * rayleigh

transmission = np.exp(-k_lambda * airmass)

# qe curve
qe_wavelength_nm = np.array([350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000])

qe_percent = np.array([20, 39, 55, 63, 72, 77, 85, 87, 86, 85, 80, 60, 38, 17])

qe_fraction = qe_percent / 100

qe_interp = interp1d(qe_wavelength_nm, qe_fraction, bounds_error=False, fill_value=0)

qe_curve = qe_interp(wavelengths_nm)

# blackbody comparison
plt.figure(figsize=(10,6))

for T in [4000, 6000, 10000]:

    spectrum = planck(wavelengths_m, T)

    spectrum /= np.max(spectrum)

    plt.plot(wavelengths_nm, spectrum, label=f"{T} K")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Flux")
plt.title("Blackbody Spectra")
plt.legend()
plt.grid(True)
plt.show()

# load filter
g_filter = filters.load_filter("lsst2023-g")

filter_interp = interp1d(g_filter.wavelength / 10, g_filter.response, bounds_error=False, fill_value=0)

filter_response = filter_interp(wavelengths_nm)

# instrument response plot
plt.figure(figsize=(10,6))

plt.plot(wavelengths_nm, filter_response, label="LSST g Filter")
plt.plot(wavelengths_nm, qe_curve, label="Detector QE")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Response")
plt.title("Instrument Response")
plt.legend()
plt.grid(True)
plt.show()

# pipeline spectrum
T = 6000

spectrum = planck(wavelengths_m, T)

spectrum /= np.trapezoid(spectrum, wavelengths_m)

atmosphere_spectrum = spectrum * transmission

filtered_spectrum = atmosphere_spectrum * filter_response

detected_spectrum = filtered_spectrum * qe_curve

plt.figure(figsize=(10,6))

plt.plot(wavelengths_nm, spectrum / np.max(spectrum), label="Blackbody")
plt.plot(wavelengths_nm, atmosphere_spectrum / np.max(spectrum), label="After Atmosphere")
plt.plot(wavelengths_nm, filtered_spectrum / np.max(spectrum), label="After Filter")
plt.plot(wavelengths_nm, detected_spectrum / np.max(spectrum), label="After QE")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative Flux")
plt.title("Pipeline Spectrum")
plt.legend()
plt.grid(True)
plt.show()

# pipeline function
def run_pipeline(T, filter_name, atmosphere=True):

    spectrum = planck(wavelengths_m, T)

    spectrum /= np.trapezoid(spectrum, wavelengths_m)

    if atmosphere:

        spectrum = spectrum * transmission

    current_filter = filters.load_filter(filter_name)

    filter_interp = interp1d(current_filter.wavelength / 10, current_filter.response, bounds_error=False, fill_value=0)

    filter_response = filter_interp(wavelengths_nm)

    spectrum *= filter_response

    spectrum *= qe_curve

    flux = np.trapezoid(spectrum, wavelengths_m)

    return flux

# multi filter comparison
filter_names = ["lsst2023-u", "lsst2023-g", "lsst2023-r", "lsst2023-i", "lsst2023-z"]

cool_fluxes = []
hot_fluxes = []

for filt in filter_names:

    cool_fluxes.append(run_pipeline(4000, filt))
    hot_fluxes.append(run_pipeline(10000, filt))

x = np.arange(len(filter_names))

width = 0.35

plt.figure(figsize=(10,6))

plt.bar(x - width/2, cool_fluxes, width, label="4000 K")
plt.bar(x + width/2, hot_fluxes, width, label="10000 K")

plt.xticks(x, ["u", "g", "r", "i", "z"])

plt.ylabel("Detected Flux")
plt.title("Filter Comparison")
plt.legend()
plt.grid(True)
plt.show()

# expected photon counts
starting_photons = 100000
exposure_time = 10

print()
print("expected photon counts")

for filt, cool_flux, hot_flux in zip(filter_names, cool_fluxes, hot_fluxes):

    cool_counts = cool_flux * starting_photons * exposure_time
    hot_counts = hot_flux * starting_photons * exposure_time

    print()
    print(f"Filter: {filt}")
    print(f"4000 K Counts = {cool_counts:.0f}")
    print(f"10000 K Counts = {hot_counts:.0f}")
    print(f"Flux Ratio = {hot_flux/cool_flux:.2f}")

# noise test
reference_flux = run_pipeline(6000, "lsst2023-g", atmosphere=True)

expected_counts = reference_flux * starting_photons * exposure_time

n_exposures = 500

photon_counts = np.random.poisson(expected_counts, n_exposures)

mean_counts = np.mean(photon_counts)
variance_counts = np.var(photon_counts)
std_counts = np.std(photon_counts)

print("noise validation")

print(f"Expected Counts = {expected_counts:.0f}")
print(f"Mean Counts = {mean_counts:.2f}")
print(f"Variance = {variance_counts:.2f}")
print(f"Standard Deviation = {std_counts:.2f}")
print(f"Variance / Mean = {variance_counts/mean_counts:.3f}")

# photon distribution
plt.figure(figsize=(8,5))
plt.hist(photon_counts, bins=25)
plt.axvline(expected_counts, color="red", linestyle="--", label="Expected Counts")

plt.xlabel("Detected Counts")
plt.ylabel("Frequency")
plt.title("Photon Count Distribution")
plt.legend()
plt.grid(True)
plt.show()

# light curve
time = np.arange(n_exposures)

plt.figure(figsize=(10,6))
plt.plot(time, photon_counts, marker=".", linestyle="-")

plt.xlabel("Exposure Number")
plt.ylabel("Detected Counts")
plt.title("Simulated Light Curve")
plt.grid(True)
plt.show()

# atmospheric loss vs temperature
temps = np.arange(3000, 12001, 1000)

losses = []

for T in temps:

    no_atm = run_pipeline(T, "lsst2023-g", atmosphere=False)
    atm = run_pipeline(T, "lsst2023-g", atmosphere=True)

    loss = 100 * (no_atm - atm) / no_atm

    losses.append(loss)

plt.figure(figsize=(10,6))
plt.plot(temps, losses, marker="o")

plt.xlabel("Temperature (K)")
plt.ylabel("Flux Loss (%)")
plt.title("Atmospheric Loss vs Temperature")
plt.grid(True)
plt.show()

# atmospheric loss for all filters
plt.figure(figsize=(10,6))
for filt in filter_names:
    losses = []

    for T in temps:
        no_atm = run_pipeline(T, filt, atmosphere=False)
        atm = run_pipeline(T, filt, atmosphere=True)

        loss = 100 * (no_atm - atm) / no_atm

        losses.append(loss)

    plt.plot(temps, losses, marker="o", label=filt[-1])

plt.xlabel("Temperature (K)")
plt.ylabel("Flux Loss (%)")
plt.title("Atmospheric Loss vs Temperature for Different Filters")
plt.legend()
plt.grid(True)
plt.show()

# differential photometry
print("differential photometry")

target_temperature = 6000
reference_temperature = 5500

filter_name = "lsst2023-g"

target_flux = run_pipeline(target_temperature, filter_name, atmosphere=True)
reference_flux = run_pipeline(reference_temperature, filter_name, atmosphere=True)

target_expected_counts = target_flux * starting_photons * exposure_time
reference_expected_counts = reference_flux * starting_photons * exposure_time

print(f"Target Temperature = {target_temperature} K")
print(f"Reference Temperature = {reference_temperature} K")
print(f"Target Flux = {target_flux:.5f}")
print(f"Reference Flux = {reference_flux:.5f}")
print(f"Target Expected Counts = {target_expected_counts:.0f}")
print(f"Reference Expected Counts = {reference_expected_counts:.0f}")
