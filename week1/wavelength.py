
import numpy as np
import matplotlib.pyplot as plt

wavelength = np.linspace(300e-9, 1100e-9, 1000)

print("Minimum wavelength:", wavelength.min())
print("Maximum wavelength:", wavelength.max())
print("Number of points:", len(wavelength))

#show plot
plt.plot(wavelength)
plt.show()
