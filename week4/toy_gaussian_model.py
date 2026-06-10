import numpy as np
import matplotlib.pyplot as plt

wavelengths_nm = np.linspace(100, 3000, 5000)

#toy gaussian model
def gaussian_filter(wavelengths, centre, width):

    response = np.exp(
        -0.5 * (
            (wavelengths - centre)/width
        )**2
    )

    return response

#create toy filter
blue_filter = gaussian_filter(
    wavelengths_nm,
    450,
    50
)

green_filter = gaussian_filter(
    wavelengths_nm,
    550,
    60
)

red_filter = gaussian_filter(
    wavelengths_nm,
    700,
    80
)

plt.figure(figsize=(10,6))

plt.plot(wavelengths_nm, blue_filter, label="Blue Filter")

plt.plot(wavelengths_nm, green_filter, label="Green Filter")

plt.plot(wavelengths_nm, red_filter, label="Red Filter")

plt.xlabel("Wavelength (nm)")

plt.ylabel("Transmission")

plt.title("Simplified Photometric Filter Responses")

plt.legend()

plt.grid(True)

plt.show()
