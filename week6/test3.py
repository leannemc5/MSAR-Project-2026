
# testing that 2 independent observations act the same
import numpy as np
import matplotlib.pyplot as plt

# input parameters
expected_counts = 70000
n_exposures = 500

lc1 = np.random.poisson(expected_counts, n_exposures)
lc2 = np.random.poisson(expected_counts, n_exposures)

# results

print("Light Curve 1")
print(f"Mean = {np.mean(lc1):.2f}")
print(f"Std = {np.std(lc1):.2f}")

print("Light Curve 2")
print(f"Mean = {np.mean(lc2):.2f}")
print(f"Std = {np.std(lc2):.2f}")

# plot
time = np.arange(n_exposures)

plt.figure(figsize=(10,6))
plt.plot(time, lc1, label="Light Curve 1", alpha=0.7)
plt.plot(time, lc2, label="Light Curve 2", alpha=0.7)

plt.xlabel("Exposure Number")
plt.ylabel("Counts")
plt.title("Repeatability Test")
plt.legend()
plt.grid(True)
plt.show()
