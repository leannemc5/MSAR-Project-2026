import numpy as np
import matplotlib.pyplot as plt
#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

wavelength = np.linspace(300e-9, 1100e-9, 1000)

# planks function
def planck(wavelength, temperature):
  num = 2*h*c**2
  exponent = (h*c)/ (wavelength*k*temperature)
  denom = wavelength**5*(np.exp(exponent) - 1)
  intensity = num/denom
  return intensity

#generate spectra 
sun = planck(wavelength, 5800)
red_star = planck(wavelength, 3500)
blue_star = planck(wavelength, 10000)

plt.plot(wavelength*1e9, red_star, label="3500 K")
plt.plot(wavelength*1e9, sun, label="5800 K")
plt.plot(wavelength*1e9, blue_star, label="10000 K")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.legend()

plt.show()
