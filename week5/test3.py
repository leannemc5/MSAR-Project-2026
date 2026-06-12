#needed libraires
import numpy as np

#brightness vs noise
count_levels = [
    1000,
    10000,
    100000
]

for counts in count_levels:

    observations = np.random.poisson(
        counts,
        1000
    )

    mean_counts = np.mean(
        observations
    )

    std_counts = np.std(
        observations
    )

    relative_noise = (
        std_counts /
        mean_counts
    )

    print()
    print(
        f"Expected Counts = {counts}"
    )
    print(
        f"Mean = {mean_counts:.2f}"
    )
    print(
        f"Std = {std_counts:.2f}"
    )
    print(
        f"Relative Noise = {relative_noise:.5f}"
    )
