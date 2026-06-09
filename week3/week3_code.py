# constants
h = 6.626e-34
c = 3e8

lambda_eff = np.array([400e-9, 550e-9, 700e-9])

#calculate photon energy
photon_energy = (h*c)/ lambda_eff
print(photon_energy)
