# import needed libraries
import numpy as np

#constants
h = 6.626e-34
c = 3.0e8
k = 1.381e-23

# planks function
def planck(wavelength, temperature):
  num = 2*h*c**2
  exponent = (h*c)/ (wavelength*k*temperature)
  denom = wavelength**5*(np.exp(exponent) - 1)
  intensity = num/denom
  return intensity
