exposure_times = [1, 5, 10, 20, 50, 100]

counts = []

for t in exposure_times:

    photon_counts = (
        flux * t
    ) / photon_energy

    counts.append(photon_counts)

import matplotlib.pyplot as plt

plt.plot(exposure_times, counts)

plt.xlabel("Exposure Time (s)")
plt.ylabel("Photon Counts")

plt.title("Photon Counts vs Exposure Time")

plt.grid(True)

plt.show()
