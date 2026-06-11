#needed libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from speclite import filters

#constants
h = 6.626e-34
c = 3e8
k = 1.381e-23

#planck Function
def planck(wavelength,T):

    exponent = h*c/(wavelength*k*T)

    return (
        2*h*c**2
    ) / (
        wavelength**5 *
        (np.exp(exponent)-1)
    )

wavelengths = np.linspace(
    300e-9,
    1000e-9,
    5000
)

wavelengths_angstrom = (
    wavelengths * 1e10
)

#blackbody
spectrum = planck(
    wavelengths,
    6000
)

#load filter
V = filters.load_filter(
    "lsst2023-g"
)

#interpolate filter
interp_filter = interp1d(
    V.wavelength,
    V.response,
    bounds_error=False,
    fill_value=0
)

filter_response = interp_filter(
    wavelengths_angstrom
)

#apply filter
filtered_spectrum = (
    spectrum *
    filter_response
)

#plot
plt.figure(figsize=(10,6))

plt.plot(
    wavelengths_angstrom,
    spectrum/np.max(spectrum),
    label="Original"
)

plt.plot(
    wavelengths_angstrom,
    filtered_spectrum/
    np.max(filtered_spectrum),
    label="Filtered"
)

plt.xlabel("Wavelength (Angstrom)")
plt.ylabel("Normalised Intensity")
plt.legend()
plt.grid()
plt.show()
