
# check to see if star stays constant over time

import numpy as np
import matplotlib.pyplot as plt

# input parameters
expected_counts = 70000
n_exposures = 500

# light curve
light_curve = np.random.poisson(expected_counts, n_exposures)

# running mean
running_mean = []
for i in range(20, n_exposures):
    running_mean.append(np.mean(light_curve[:i]))

# results
print(f"Expected Counts = {expected_counts}")
print(f"Final Running Mean = {running_mean[-1]:.2f}")

# plot
plt.figure(figsize=(10,6))
plt.plot(np.arange(20, n_exposures), running_mean, label="Running Mean")
plt.axhline(expected_counts, color="red", linestyle="--", label="Expected Counts")

plt.xlabel("Exposure Number")
plt.ylabel("Mean Counts")
plt.title("Running Mean Stability")
plt.legend()
plt.grid(True)
plt.show()
