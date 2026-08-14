import numpy as np
import matplotlib.pyplot as plt

# scatter comparison
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

target_std = np.std(target_curve)
reference_std = np.std(reference_curve)
ratio_std = np.std(ratio_curve)

print(f"Target Std = {target_std:.2f}")
print(f"Reference Std = {reference_std:.2f}")
print(f"Ratio Std = {ratio_std:.5f}")

plt.figure(figsize=(8,5))
plt.bar(["Target", "Reference", "Ratio"], [target_std, reference_std, ratio_std])

plt.ylabel("Standard Deviation")
plt.title("Scatter Comparison")
plt.grid(True)
plt.show()
