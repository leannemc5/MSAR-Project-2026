import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import h, c, k
from scipy.interpolate import interp1d
from speclite import filters

# plancks function
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


# wavelength grid

wavelengths_nm = np.linspace(
    300,
    1100,
    5000
)

wavelengths_m = wavelengths_nm * 1e-9

# stellar spectrum

T = 6000

spectrum = planck(
    wavelengths_m,
    T
)

# Normalise spectrum area

spectrum = (
    spectrum
    /
    np.trapezoid(
        spectrum,
        wavelengths_m
    )
)

# filter response
g_filter = filters.load_filter(
    "lsst2023-g"
)

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
    spectrum
    *
    filter_response
)

# detector qe
qe_wavelength_nm = np.array([
    350,400,450,500,550,
    600,650,700,750,800,
    850,900,950,1000
])

qe_percent = np.array([
    20,39,55,63,72,
    77,85,87,86,85,
    80,60,38,17
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
    filtered_spectrum
    *
    qe_curve
)

#integrate flux
total_flux = np.trapezoid(
    detected_spectrum,
    wavelengths_m
)

print(
    f"Integrated Flux = {total_flux:.4f}"
)

#photon counts
starting_photons = 100000

exposure_time = 10

expected_counts = (
    total_flux
    *
    starting_photons
    *
    exposure_time
)

print(f"Expected Photon Counts = {expected_counts:.0f}")

# poisson and read noise
n_exposures = 500

poisson_counts = np.random.poisson(
    expected_counts,
    n_exposures
)

read_noise_sigma = 100

read_noise = np.random.normal(
    0,
    read_noise_sigma,
    n_exposures
)

light_curve = (poisson_counts + read_noise)

# stats
mean_counts = np.mean(
    light_curve
)

variance_counts = np.var(
    light_curve
)

std_counts = np.std(
    light_curve
)

print(f"Mean Counts = {mean_counts:.2f}")
print(f"Variance = {variance_counts:.2f}")
print(f"Standard Deviation = {std_counts:.2f}")
print(f"Variance / Mean = {variance_counts/mean_counts:.3f}")


time = np.arange(
    n_exposures
)

# pipeline plot

plt.figure(figsize=(10,6))

plt.plot(
    wavelengths_nm,
    spectrum / np.max(spectrum),
    label="Blackbody"
)

plt.plot(
    wavelengths_nm,
    filtered_spectrum / np.max(spectrum),
    label="After Filter"
)

plt.plot(
    wavelengths_nm,
    detected_spectrum / np.max(spectrum),
    label="After QE"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Relative Flux")

plt.title(
    f"Pipeline Spectrum ({T} K)"
)

plt.legend()

plt.grid(True)

plt.show()


# QE curve
plt.figure(figsize=(8,5))

plt.plot(
    qe_wavelength_nm,
    qe_percent,
    marker="o"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Quantum Efficiency (%)")

plt.title(
    "Detector QE"
)

plt.grid(True)

plt.show()


#photon dist
plt.figure(figsize=(8,5))

plt.hist(
    light_curve,
    bins=25
)

plt.axvline(
    expected_counts,
    color="red",
    linestyle="--",
    label="Expected Counts"
)

plt.xlabel("Detected Counts")
plt.ylabel("Frequency")

plt.title(
    "Photon Count Distribution"
)

plt.legend()

plt.grid(True)

plt.show()


# week 6 light curve
plt.figure(figsize=(10,6))

plt.plot(
    time,
    light_curve,
    marker=".",
    linestyle="-"
)

plt.xlabel("Exposure Number")
plt.ylabel("Detected Counts")

plt.title(
    "Simulated Light Curve"
)

plt.grid(True)

plt.show()

#noise over time

running_std = []

for i in range(20, n_exposures):

    running_std.append(
        np.std(
            light_curve[:i]
        )
    )

plt.figure(figsize=(10,6))

plt.plot(
    np.arange(
        20,
        n_exposures
    ),
    running_std
)

plt.xlabel("Exposure Number")
plt.ylabel("Running Standard Deviation")

plt.title(
    "Noise Stability Over Time"
)

plt.grid(True)

plt.show()
