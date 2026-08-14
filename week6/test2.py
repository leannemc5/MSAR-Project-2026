
# test to see what exposure times give more stable estimates

import numpy as np
import matplotlib.pyplot as plt

# input parameters

expected_counts = 70000
exposure_numbers = [50, 200, 1000]

means = []
stds = []

for n in exposure_numbers:

    light_curve = np.random.poisson(expected_counts, n)

    means.append(np.mean(light_curve))
    stds.append(np.std(light_curve))

    print()
    print(f"Exposures = {n}")
    print(f"Mean = {means[-1]:.2f}")
    print(f"Std = {stds[-1]:.2f}")

# plot

plt.figure(figsize=(8,5))

plt.plot(exposure_numbers, means, marker="o")

plt.axhline(expected_counts, color="red", linestyle="--", label="Expected Counts")

plt.xlabel("Number of Exposures")
plt.ylabel("Mean Counts")
plt.title("Mean Brightness vs Observation Length")
plt.legend()
plt.grid(True)
plt.show()
