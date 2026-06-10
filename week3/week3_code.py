#needed libraries
import numpy as np
import matplotlib.pyplot as plt

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

#plancks function
def planck(wavelength, temperature):

    exponent = (h * c) / (wavelength * k * temperature)

    intensity = (
        (2 * h * c**2)
        /
        (wavelength**5 * (np.exp(exponent) - 1))
    )

    return intensity

wavelengths = np.linspace(100e-9, 3000e-9, 5000)

#stellar temps
temperature = 6000

spectrum = planck(wavelengths, temperature)

#total flux
total_flux = np.trapezoid(spectrum, wavelengths)

print("Integrated Flux =", total_flux)

#wavelength
lambda_eff = 550e-9

#photon energy
photon_energy = (h * c) / lambda_eff

print("Photon Energy =", photon_energy)

#exposure time
exposure_time = 100

#photon counts
photon_counts = (
    total_flux * exposure_time
) / photon_energy

print("Photon Counts =", photon_counts)
