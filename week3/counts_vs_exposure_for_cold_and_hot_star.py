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

    numerator = 2 * h * c**2

    denominator = wavelength**5 * (np.exp(exponent) - 1)

    intensity = numerator / denominator

    return intensity


wavelengths = np.linspace(100e-9, 3000e-9, 5000)

# convert to nm for plotting
wavelengths_nm = wavelengths * 1e9
lambda_eff = 550e-9

# photon energy
photon_energy = (h * c) / lambda_eff

#stellar temps
temperatures = [3000, 10000]

# colours for plotting
colours = {
    3000: "red",
    10000: "blue"
}

exposure_times = np.linspace(1, 100, 50)

plt.figure(figsize=(10,6))

for T in temperatures:
    spectrum = planck(wavelengths, T)

    #integrate total flux
    total_flux = np.trapezoid(
        wavelengths
    )

    #calculates count for each exposure
    counts = []

    for t in exposure_times:

        photon_counts = (
            total_flux * t
        ) / photon_energy

        counts.append(photon_counts)

    #plot counts vs exposure time
    plt.plot(
        exposure_times,
        counts,
        color=colours[T],
        label=f"{T} K Star"
    )
plt.xlabel("Exposure Time (s)")
plt.ylabel("Photon Counts")
plt.title("Photon Counts vs Exposure Time")
plt.legend()
plt.grid(True)
plt.show()
