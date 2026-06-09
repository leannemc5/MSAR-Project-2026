#needed libraries
import numpy as np
import matplotlib.pyplot as plt

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

# Wien constant
b = 2.9e-3   # m K

def planck(wavelength, temperature):

    exponent = (h * c) / (wavelength * k * temperature)

    numerator = 2 * h * c**2

    denominator = wavelength**5 * (np.exp(exponent) - 1)

    return numerator / denominator

wavelengths = np.linspace(100e-9, 3000e-9, 5000)
temperatures = [3000, 5000, 7000, 10000]

#plot
plt.figure(figsize=(10,6))
for T in temperatures:

    spectrum = planck(wavelengths, T)

    # normalise for plotting
    spectrum = spectrum / np.max(spectrum)

    # Find numerical peak
    peak_index = np.argmax(spectrum)
    peak_wavelength = wavelengths[peak_index]

    # wein
    predicted_peak = b / T
    print(f"Temperature = {T} K")

    print(
        f"Numerical Peak = {peak_wavelength*1e9:.1f} nm"
    )

    print(
        f"Wien Prediction = {predicted_peak*1e9:.1f} nm"
    )

    print()

    plt.plot(
        wavelengths * 1e9,
        spectrum,
        label=f"{T} K"
    )


    plt.axvline(
        peak_wavelength * 1e9,
        linestyle="--"
    )

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Intensity")
plt.title("Blackbody Spectra and Wien's Law Peaks")
plt.legend()
plt.grid(True)
plt.show()
