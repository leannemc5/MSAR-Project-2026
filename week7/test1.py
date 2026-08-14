import numpy as np
import matplotlib.pyplot as plt

# test of atmospheric extinction with airmass

airmass = np.linspace(1, 3, 100)

k = 0.2

transmission = np.exp(-k * airmass)

plt.figure(figsize=(8,5))

plt.plot(airmass, transmission)

plt.xlabel("Airmass")
plt.ylabel("Transmission")
plt.title("Atmospheric Extinction")
plt.grid(True)
plt.show()
