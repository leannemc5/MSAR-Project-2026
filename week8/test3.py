import numpy as np
import matplotlib.pyplot as plt

# different brightness target and reference
expected_counts_target = 60000
expected_counts_reference = 40000

n_exposures = 500

cloud_factor = np.random.normal(0.90, 0.03, n_exposures)

cloud_factor = np.clip(cloud_factor, 0, 1)

target_curve = []
reference_curve = []
ratio_curve = []

for transparency in cloud_factor:

    target = np.random.poisson(expected_counts_target * transparency)
    reference = np.random.poisson(expected_counts_reference * transparency)

    target_curve.append(target)
    reference_curve.append(reference)
    ratio_curve.append(target/reference)

target_curve = np.array(target_curve)
reference_curve = np.array(reference_curve)
ratio_curve = np.array(ratio_curve)

print(f"Expected Ratio = {60000/40000:.3f}")
print(f"Measured Ratio = {np.mean(ratio_curve):.3f}")

plt.figure(figsize=(10,6))
plt.hist(ratio_curve, bins=30)

plt.xlabel("Target / Reference")
plt.ylabel("Frequency")
plt.title("Distribution of Differential Ratio")
plt.grid(True)
plt.show()

plt.figure(figsize=(10,6))
plt.plot(ratio_curve)
plt.xlabel("Exposure")
plt.ylabel("Target / Reference")
plt.title("Differential Light Curve (Different Brightness)")
plt.grid(True)
plt.show()
