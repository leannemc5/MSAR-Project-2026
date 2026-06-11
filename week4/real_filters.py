import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from speclite import filters

#constants
h = 6.626e-34
c = 3e8
k = 1.381e-23

#planck function
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

#normalise the spectrum
spectrum = (
    spectrum /
    np.trapezoid(
        spectrum,
        wavelengths_m
    )
)

print(
    "Spectrum area:",
    np.trapezoid(
        spectrum,
        wavelengths_m
    )
)

#real filter
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

#add filter
filtered_spectrum = (
    spectrum *
    filter_response
)

print(
    "Filtered area:",
    np.trapezoid(
        filtered_spectrum,
        wavelengths_m
    )
)

#detector qe data
qe_wavelength_nm = np.array([
    350, 400, 450, 500, 550, 600,
    650, 700, 750, 800, 850, 900,
    950, 1000
])

qe_percent = np.array([
    20, 39, 55, 63, 72, 77,
    85, 87, 86, 85, 80, 60,
    38, 17
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

#add qe
detected_spectrum = (
    filtered_spectrum *
    qe_curve
)

print(
    "Detected area:",
    np.trapezoid(
        detected_spectrum,
        wavelengths_m
    )
)

#integrate flux
total_flux = np.trapezoid(
    detected_spectrum,
    wavelengths_m
)

print(
    f"\nIntegrated Flux = {total_flux:.4f}"
)


# toy photon count model
starting_photons = 100000

exposure_time = 10

photon_counts = (
    total_flux *
    starting_photons *
    exposure_time
)

print(
    f"Detected Photons = {photon_counts:.0f}"
)

#plot spectrum
plt.figure(figsize=(10,6))

plt.plot(
    wavelengths_nm,
    spectrum / np.max(spectrum),
    label="Original Spectrum"
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
plt.ylabel("Relative Intensity")
plt.title("Spectrum Through Filter and Detector")
plt.legend()
plt.grid(True)
plt.show()
