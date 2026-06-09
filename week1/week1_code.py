#needed libraries
import numpy as np
import matplotlib.pyplot as plt

#constants
h = 6.626e-34
c = 2.998e8
k = 1.381e-23

#sample range of wavelengths
wavelength = np.linspace(100e-9, 3000e-9, 5000)

# planks function
def planck(wavelength, temperature):
  num = 2*h*c**2
  exponent = (h*c)/ (wavelength*k*temperature)
  denom = wavelength**5*(np.exp(exponent) - 1)
  intensity = num/denom
  return intensity

#stellar temps
temperatures = [3000, 5000, 7000, 10000]

#plot
plt.figure(figsize=(10,6))

for T in temperatures:
  spectrum = planck(wavelength, T)

  # normalise spectrafor easier comparison
  spectrum = spectrum/np.max(spectrum)

  plt.plot(
      wavelength* 1e9,
      spectrum,
      label=f"{T} K"
  )

plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalised Intensity")
plt.title("Blackbody Spectra")
plt.legend()
plt.grid()
plt.show()
