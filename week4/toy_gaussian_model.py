#needed libraires
import numpy as np
import matplotlib.pyplot as plt

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

#plancks function
def planck(wavelength, temperature):

    exponent = (h*c)/(wavelength*k*temperature)

    numerator = 2*h*c**2

    denominator = (
        wavelength**5 *
        (np.exp(exponent)-1)
    )

    intensity = numerator / denominator

    return intensity

#toy gaussian
def gaussian_filter(wavelengths, centre, width):

    response = np.exp(
        -0.5 * (
            (wavelengths - centre)/width
        )**2
    )

    return response

wavelengths = np.linspace(
    100e-9,
    3000e-9,
    5000
)

wavelengths_nm = wavelengths * 1e9

#spectra
cool_star = planck(wavelengths, 3000)

hot_star = planck(wavelengths, 10000)

#normalised
cool_star = cool_star / np.max(cool_star)

hot_star = hot_star / np.max(hot_star)

#toy filters
blue_filter = gaussian_filter(
    wavelengths_nm,
    450,
    50
)

red_filter = gaussian_filter(
    wavelengths_nm,
    700,
    80
)

#add filters to spectra
cool_blue = cool_star * blue_filter
cool_red = cool_star * red_filter
hot_blue = hot_star * blue_filter
hot_red = hot_star * red_filter

#plot
plt.figure(figsize=(12,7))

plt.plot(
    wavelengths_nm,
    cool_star,
    label="Cool Star (3000 K)"
)

plt.plot(
    wavelengths_nm,
    hot_star,
    label="Hot Star (10000 K)"
)

plt.plot(
    wavelengths_nm,
    blue_filter,
    '--',
    label="Blue Filter"
)

plt.plot(
    wavelengths_nm,
    red_filter,
    '--',
    label="Red Filter"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalized Intensity")
plt.title("Stellar Spectra with Filter Responses")
plt.legend()
plt.grid(True)
plt.show()
