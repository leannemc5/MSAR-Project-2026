#needed libraries
import numpy as np
from scipy.interpolate import interp1d

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

#plancks function
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

def calculate_photon_counts(
    temperature,
    exposure_time
):

    wavelengths_nm = np.linspace(
        300,
        1000,
        5000
    )

    wavelengths_m = (
        wavelengths_nm * 1e-9
    )

    #make stellar spectrum
    spectrum = planck(
        wavelengths_m,
        temperature
    )

    #add detector QE
    qe_interp = interp1d(
        qe_wavelength_nm,
        qe,
        bounds_error=False,
        fill_value=0
    )

    qe_curve = qe_interp(
        wavelengths_nm
    )

    detected_spectrum = (
        spectrum *
        qe_curve
    )

    #integrate flux
    total_flux = np.trapezoid(
        detected_spectrum,
        wavelengths_m
    )

    #reference wavelength
    wavelength_ref = 550e-9

    photon_energy = (
        h * c
    ) / wavelength_ref

    photon_rate = (
        total_flux /
        photon_energy
    )

    photon_counts = (
        photon_rate *
        exposure_time
    )

    return photon_counts

if __name__ == "__main__":

    counts = calculate_photon_counts(
        temperature=6000,
        exposure_time=10
    )

    print(
        f"Photon Counts: {counts:.3e}"
    )
