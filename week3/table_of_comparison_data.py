#define list of values to test
fluxes = [1e-9, 2e-9]
exposures = [5, 10, 20]
wavelengths = [400e-9, 700e-9]

#table header
print(f"{'Flux (W)':<10} | {'Exp (s)':<7} | {'Wavelength':<12} | {'Photon Counts'}")
print("-" * 55)

#loop through every single combo
for flux in fluxes:
    for exp in exposures:
        for wl in wavelengths:

            # calculate photon energy
            photon_energy = (h * c) / wl

            #calculate photon counts
            photon_counts = (flux * exp) / photon_energy

            #print results
            print(f"{flux:<10.1e} | {exp:<7} | {wl*1e9:<8.0f} nm | {photon_counts:,.2f}")
