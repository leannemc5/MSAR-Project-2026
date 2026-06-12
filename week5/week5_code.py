#needed libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from speclite import filters

#constants
h = 6.626e-34
c = 3e8
k = 1.381e-23

#plankcs function
def planck(wavelength, T):

    exponent = (h * c) / (wavelength * k * T)

    intensity = (
        (2 * h * c**2)
        /
        (wavelength**5)
    ) / (
        np.exp(exponent) - 1
    )

    return intensity

wavelengths_nm = np.linspace(
    300,
    1100,
    5000
)

wavelengths_m = wavelengths_nm * 1e-9

T = 6000

spectrum = planck(
    wavelengths_m,
    T
)

#mormalise spectrum
spectrum = (
    spectrum
    /
    np.trapezoid(
        spectrum,
        wavelengths_m
    )
)

#filter response
g_filter = filters.load_filter(
    "lsst2023-g"
)

#filter interpolate
filter_interp = interp1d(
    g_filter.wavelength / 10,
    g_filter.response,
    bounds_error=False,
    fill_value=0
)

filter_response = filter_interp(
    wavelengths_nm
)

filtered_spectrum = (
    spectrum *
    filter_response
)

#detector qe data
qe_wavelength_nm = np.array([
    350, 400, 450, 500, 550,
    600, 650, 700, 750, 800,
    850, 900, 950, 1000
])

qe_percent = np.array([
    20, 39, 55, 63, 72,
    77, 85, 87, 86, 85,
    80, 60, 38, 17
])

qe_fraction = qe_percent / 100

qe_interp = interp1d(
    qe_wavelength_nm,
    qe_fraction,
    bounds_error=False,
    fill_value=0
)

qe_curve = qe_interp(
    wavelengths_nm
)

detected_spectrum = (
    filtered_spectrum *
    qe_curve
)

#integarate flux
total_flux = np.trapezoid(
    detected_spectrum,
    wavelengths_m
)

print()
print(
    f"Integrated Flux = {total_flux:.4f}"
)

#photon counts
starting_photons = 1000000

exposure_time = 10

expected_counts = (
    total_flux
    *
    starting_photons
    *
    exposure_time
)

print(
    f"Expected Photon Counts = {expected_counts:.0f}"
)

#noise model
n_measurements = 500

noisy_counts = np.random.poisson(
    expected_counts,
    n_measurements
)

mean_counts = np.mean(
    noisy_counts
)

variance_counts = np.var(
    noisy_counts
)

print()
print(
    f"Mean Counts = {mean_counts:.2f}"
)

print(
    f"Variance = {variance_counts:.2f}"
)

print(
    f"Variance / Mean = {variance_counts/mean_counts:.3f}"
)


#spectrum plot
plt.figure(figsize=(10,6))

plt.plot(
    wavelengths_nm,
    spectrum / np.max(spectrum),
    label="Blackbody"
)

plt.plot(
    wavelengths_nm,
    filtered_spectrum / np.max(filtered_spectrum),
    label="After Filter"
)

plt.plot(
    wavelengths_nm,
    detected_spectrum / np.max(detected_spectrum),
    label="After QE"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Intensity")

plt.title(f"Pipeline Spectrum ({T} K)"
)

plt.legend()
plt.grid(True)
plt.show()

#qe curve plot
plt.figure(figsize=(8,5))

plt.plot(
    qe_wavelength_nm,
    qe_percent,
    marker="o"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Quantum Efficiency (%)")

plt.title("Detector Quantum Efficiency")
plt.grid(True)
plt.show()

#histogram
plt.hist(
    noisy_counts,
    bins=25,
    alpha=0.7,
    label="Noisy Measurements"
)

plt.axvline(
    expected_counts,
    color="red",
    linestyle="--",
    label="Expected Counts"
)

plt.xlabel("Detected Photon Counts")
plt.ylabel("Frequency")

plt.title(
    "Poisson Noise Distribution"
)
plt.legend()
plt.grid(True)
plt.show()
