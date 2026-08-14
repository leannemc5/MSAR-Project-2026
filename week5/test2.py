#needed libraries
import numpy as np
import matplotlib.pyplot as plt

#histogram
expected_counts = 70000

observed_counts = np.random.poisson(
    expected_counts,
    1000
)

plt.figure(figsize=(8,5))

plt.hist(
    observed_counts,
    bins=25
)

plt.xlabel("Photon Counts")
plt.ylabel("Frequency")
plt.title("Photon Count Distribution")
plt.grid(True)
plt.show()
