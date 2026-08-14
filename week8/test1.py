import numpy as np
import matplotlib.pyplot as plt

#identical target and reference stars
expected_counts_target = 60000
expected_counts_reference = 60000

n_exposures = 500
read_noise_sigma = 100

cloud_factor = np.random.normal(0.90, 0.03, n_exposures)

cloud_factor = np.clip(cloud_factor, 0, 1)

target_curve = []
reference_curve = []
ratio_curve = []

for transparency in cloud_factor:

    target = np.random.poisson(expected_counts_target * transparency)
    target += np.random.normal(0, read_noise_sigma)

    reference = np.random.poisson(expected_counts_reference * transparency)
    reference += np.random.normal(0, read_noise_sigma)

    target_curve.append(target)
    reference_curve.append(reference)
    ratio_curve.append(target/reference)

target_curve = np.array(target_curve)
reference_curve = np.array(reference_curve)
ratio_curve = np.array(ratio_curve)

print(f"Target Mean = {np.mean(target_curve):.2f}")
print(f"Reference Mean = {np.mean(reference_curve):.2f}")
print(f"Ratio Mean = {np.mean(ratio_curve):.4f}")
print(f"Ratio Std = {np.std(ratio_curve):.5f}")

time = np.arange(n_exposures)

plt.figure(figsize=(10,6))
plt.plot(time, target_curve, label="Target")
plt.plot(time, reference_curve, label="Reference")
plt.xlabel("Exposure")
plt.ylabel("Counts")
plt.title("Target and Reference Light Curves")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10,6))
plt.plot(time, ratio_curve, color="black")
plt.xlabel("Exposure")
plt.ylabel("Target / Reference")
plt.title("Differential Light Curve")
plt.grid(True)
plt.show()
