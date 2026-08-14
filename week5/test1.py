#needed libraries
import numpy as np

#variance vs mean
expected_counts = 70000

n_observations = 1000

observed_counts = np.random.poisson(
    expected_counts,
    n_observations
)

mean_counts = np.mean(
    observed_counts
)

variance_counts = np.var(
    observed_counts
)

std_counts = np.std(
    observed_counts
)

print()
print(
    f"Mean = {mean_counts:.2f}"
)
print(
    f"Variance = {variance_counts:.2f}"
)
print(
    f"Standard Deviation = {std_counts:.2f}"
)
print()
print(
    f"Variance / Mean = {variance_counts/mean_counts:.3f}"
)
