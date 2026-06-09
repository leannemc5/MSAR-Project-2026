# constants
h = 6.626e-34
c = 3e8

#fake flux
flux = 1e-9
exposure_time = 10

lambda_eff = np.array([400e-9, 550e-9, 700e-9])

#calculate photon energy
photon_energy = (h*c)/ lambda_eff

photon_counts = (
    flux*exposure_time
)/photon_energy
print(photon_counts)
