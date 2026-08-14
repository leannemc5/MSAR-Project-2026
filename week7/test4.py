import numpy as np
import matplotlib.pyplot as plt

expected_counts = 100000

n_exposures = 500

cloud_factor = (1 + np.random.normal(0, 0.03, n_exposures))

light_curve = []

for transparency in cloud_factor:

    counts = expected_counts * transparency

    measured = np.random.poisson(counts)

    light_curve.append(measured)

plt.figure(figsize=(10,6))

plt.plot(light_curve)

plt.xlabel("Exposure")
plt.ylabel("Counts")
plt.title("Light Curve With Transparency Variations")
plt.grid(True)
plt.show()
