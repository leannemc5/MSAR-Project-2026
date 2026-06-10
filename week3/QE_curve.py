#needed libraires
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

#plancks gunction
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

T = 6000

spectrum = planck(
    wavelengths_m,
    T
)

# normalise spectrum
spectrum = spectrum / np.max(spectrum)

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

#add detector spectrum
detected_spectrum = (
    spectrum *
    qe_curve
)

# QE curve
plt.figure(figsize=(10,6))
plt.plot(
    wavelengths_nm,
    qe_curve
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Quantum Efficiency")
plt.title("Detector QE Curve")
plt.grid(True)
plt.show()

#plot og spectrum against new spectrum
plt.figure(figsize=(10,6))

plt.plot(
    wavelengths_nm,
    spectrum,
    label="Original Spectrum"
)

plt.plot(
    wavelengths_nm,
    detected_spectrum,
    label="After Detector QE"
)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Intensity")
plt.title("Effect of Detector Quantum Efficiency")
plt.legend()
plt.grid(True)
plt.show()
