#needed libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

#planck function
def planck(wavelength, T):

    exponent = (h * c) / (wavelength * k * T)

    intensity = (
        (2 * h * c**2)
        /
        (
            wavelength**5
            *
            (np.exp(exponent) - 1)
        )
    )

    return intensity

wavelengths_nm = np.linspace(
    300,
    1000,
    5000
)

wavelengths_m = wavelengths_nm * 1e-9

#detector de data
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

qe = qe_percent / 100

#interpolate qe
qe_interp = interp1d(
    qe_wavelength_nm,
    qe,
    bounds_error=False,
    fill_value=0
)

qe_curve = qe_interp(
    wavelengths_nm
)

#plots of QE curves
temperatures = [
    3000,
    6000,
    10000
]

plt.figure(figsize=(10,6))

for T in temperatures:

    spectrum = planck(
        wavelengths_m,
        T
    )

    #normalise blackbody
    spectrum = (
        spectrum /
        np.max(spectrum)
    )

    #add QE
    detected_spectrum = (
        spectrum *
        qe_curve
    )

    # normalise
    detected_spectrum = (
        detected_spectrum /
        np.max(detected_spectrum)
    )

    plt.plot(
        wavelengths_nm,
        detected_spectrum,
        label=f"{T} K"
    )

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Intensity")
plt.title("Detector QE Corrected Spectra")
plt.legend()
plt.grid(True)
plt.show()
