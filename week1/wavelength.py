
# import needed libraries
import numpy as np

#sample wavelengths in metres
wavelength = np.linspace(300e-9, 1000e-9, 1000)

print("Minimum wavelength:", wavelength.min())
print("Maximum wavelength:", wavelength.max())
print("Number of points:", len(wavelength))
