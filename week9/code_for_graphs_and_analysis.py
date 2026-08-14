# blackbody spectra
temperatures = [4000, 6000, 10000]

plt.figure(figsize=(10, 6))

for T in temperatures:
    spectrum = planck(wavelengths_m, T)
    spectrum_norm = spectrum / np.max(spectrum)
    peak = wien_peak(T)

    plt.plot(wavelengths_nm, spectrum_norm, linewidth=2, label=f'{T} K')
    plt.axvline(x=peak, linestyle='--', alpha=0.5)

plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalised Flux')
plt.title('Blackbody Spectra at Different Temperatures')
plt.legend()
plt.grid(True)
plt.xlim([300, 1100])
plt.show()

# atmospheric transmission

airmasses = [1.0, 1.5, 2.0, 2.5]

plt.figure(figsize=(10, 6))
for am in airmasses:
    transmission = atmosphere_transmission(wavelengths_nm, am)
    plt.plot(wavelengths_nm, transmission, linewidth=2, label=f'airmass = {am}')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Atmospheric Transmission')
plt.title('Atmospheric Transmission vs Airmass')
plt.legend()
plt.grid(True)
plt.show()

# noise distribution

flux_test = run_pipeline(6000, "lsst2023-g", True)
expected_counts = flux_test * starting_photons * exposure_time

n_samples = 1000
clean_counts = np.random.poisson(expected_counts, n_samples)
noisy_counts = clean_counts + np.random.normal(0, read_noise_sigma, n_samples)

plt.figure(figsize=(10, 6))

plt.hist(noisy_counts, bins=30, alpha=0.7, edgecolor='black')
plt.axvline(x=expected_counts, linestyle='--', linewidth=2, label=f'Expected: {expected_counts:.0f}')
plt.axvline(x=np.mean(noisy_counts), color='r', linestyle='--', linewidth=2, label=f'Mean: {np.mean(noisy_counts):.0f}')

plt.xlabel('Detected Counts')
plt.ylabel('Frequency')
plt.title('Photon Count Distribution')
plt.legend()
plt.grid(True)
plt.show()

print(f"Expected counts: {expected_counts:.0f}")
print(f"Mean detected: {np.mean(noisy_counts):.1f}")
print(f"Std detected: {np.std(noisy_counts):.1f}")
print(f"Variance/Mean ratio: {np.var(noisy_counts)/np.mean(noisy_counts):.3f}")

# single light curve

flux_test = run_pipeline(6000, "lsst2023-g", True)
expected_counts = flux_test * starting_photons * exposure_time

light_curve_data, cloud = light_curve(expected_counts, n_exposures)
time = np.arange(n_exposures)

plt.figure(figsize=(10, 6))

plt.plot(time, light_curve_data, linewidth=0.5)
plt.axhline(y=expected_counts, color='r', linestyle='--', linewidth=2, label=f'Expected ({expected_counts:.0f})')

plt.xlabel('Exposure Number')
plt.ylabel('Detected Counts')
plt.title('Simulated Light Curve (6000 K, g-band)')
plt.legend()
plt.grid(True)
plt.show()

# diff photometry for matched temps

result = differential_photometry(6000, 6000, "lsst2023-g", n_exposures, exposure_time, starting_photons)
time = np.arange(n_exposures)

fig, axes = plt.subplots(3, 1, figsize=(10, 12))

# target and reference light curves

axes[0].plot(time, result['target_lc'], linewidth=0.5, label='Target (6000K)')
axes[0].plot(time, result['reference_lc'], linewidth=0.5, label='Reference (6000K)')
axes[0].set_xlabel('Exposure Number')
axes[0].set_ylabel('Detected Counts')
axes[0].set_title('Target and Reference Light Curves ($\Delta$T = 0 K)')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# differential light curve

axes[1].plot(time, result['diff_lc'], linewidth=0.5)
axes[1].axhline(y=1.0, linestyle='--', linewidth=2, label='Ratio = 1.0')
axes[1].axhline(y=np.mean(result['diff_lc']), linestyle='--', linewidth=1,
                label=f'Mean = {np.mean(result["diff_lc"]):.5f}')
axes[1].set_xlabel('Exposure Number')
axes[1].set_ylabel('Target / Reference')
axes[1].set_title('Differential Light Curve ($\Delta$T = 0 K)')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

# histogram

axes[2].hist(result['diff_lc'], bins=30, alpha=0.7, edgecolor='black')
axes[2].axvline(x=1.0, linestyle='--', linewidth=2, label='Expected')
axes[2].axvline(x=np.mean(result['diff_lc']), linestyle='--', linewidth=2,
                label=f'Mean = {np.mean(result["diff_lc"]):.5f}')
axes[2].set_xlabel('Target / Reference Ratio')
axes[2].set_ylabel('Frequency')
axes[2].set_title('Distribution of Differential Measurements ($\Delta$T = 0 K)')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print(f"Differential mean: {np.mean(result['diff_lc']):.5f}")
print(f"Differential std: {np.std(result['diff_lc']):.5f}")
print(f"Relative error: {np.std(result['diff_lc'])/np.mean(result['diff_lc']):.5f}")

# spectral mismatch mean ratio vs delta t

target_temp = 6000
delta_t_range = np.arange(-4000, 4001, 500)

mean_ratios = []
std_ratios = []
rel_errors = []
delta_ts = []

for delta_t in delta_t_range:
    ref_temp = target_temp - delta_t
    if ref_temp > 0:
        result = differential_photometry(target_temp, ref_temp, "lsst2023-g",
                                        n_exposures, exposure_time, starting_photons)

        delta_ts.append(delta_t)
        mean_ratios.append(np.mean(result['diff_lc']))
        std_ratios.append(np.std(result['diff_lc']))
        rel_errors.append(np.std(result['diff_lc']) / np.mean(result['diff_lc']))

plt.figure(figsize=(10, 6))

plt.plot(delta_ts, mean_ratios, 'o-', linewidth=2, markersize=6)
plt.axhline(y=1.0, linestyle='--', alpha=0.5, label='Ratio = 1.0 (no mismatch)')

plt.xlabel('Temperature Difference (Target - Reference) [K]')
plt.ylabel('Mean Flux Ratio (Target/Reference)')
plt.title('Spectral Mismatch: Mean Ratio vs $\Delta$T')
plt.legend()
plt.grid(True)
plt.show()

# spectral mismatch relative error vs delta t

abs_delta_ts = [abs(dt) for dt in delta_ts]

plt.figure(figsize=(10, 6))

plt.plot(delta_ts, rel_errors, 'o-', linewidth=2, markersize=6)
plt.axhline(y=rel_errors[abs_delta_ts.index(0)], color='r', linestyle='--', alpha=0.5,
            label=f'Baseline error ($\Delta$T = 0)')

plt.xlabel('Temperature Difference [K]')
plt.ylabel('Relative Error')
plt.title('Precision Degradation with Spectral Mismatch')
plt.legend()
plt.grid(True)
plt.show()

# filter comparison for mismatch

filters_to_plot = ["lsst2023-u", "lsst2023-g", "lsst2023-r", "lsst2023-i", "lsst2023-z"]
delta_ts_filter = [1000, 2000, 3000, 4000]

filter_errors = {}
for filt in filters_to_plot:
    errors = []
    for delta_t in delta_ts_filter:
        ref_temp = 6000 - delta_t
        result = differential_photometry(6000, ref_temp, filt,
                                        n_exposures, exposure_time, starting_photons)
        errors.append(np.std(result['diff_lc']) / np.mean(result['diff_lc']))
    filter_errors[filt] = errors

plt.figure(figsize=(10, 6))

for filt in filters_to_plot:
    plt.plot(delta_ts_filter, filter_errors[filt], 'o-', linewidth=2, markersize=6, label=filt[-1])

plt.xlabel('Temperature Difference (K)')
plt.ylabel('Relative Error')
plt.title('Spectral Mismatch Error Across Different Filters')
plt.legend()
plt.grid(True)
plt.show()

print("\nRelative error vs temperature difference for each filter:")
print(f"{'Filter':>8} {'ΔT=1000K':>12} {'ΔT=2000K':>12} {'ΔT=3000K':>12} {'ΔT=4000K':>12}")
print("-" * 56)
for filt in filters_to_plot:
    print(f"{filt[-1]:>8} {filter_errors[filt][0]:12.5f} {filter_errors[filt][1]:12.5f} "
          f"{filter_errors[filt][2]:12.5f} {filter_errors[filt][3]:12.5f}")

# temporal binning

result_bin = differential_photometry(6000, 4000, "lsst2023-g",
                                    n_exposures, exposure_time, starting_photons)

bin_sizes = [1, 2, 5, 10, 20, 50]
binned_stds = []

for b in bin_sizes:
    binned = bin_data(result_bin['diff_lc'], b)
    binned_stds.append(np.std(binned))

plt.figure(figsize=(10, 6))

plt.plot(bin_sizes, binned_stds, 'o-', linewidth=2, markersize=8)

expected_improvement = [np.std(result_bin['diff_lc'])/np.sqrt(b) for b in bin_sizes]
plt.plot(bin_sizes, expected_improvement, '--', linewidth=2, alpha=0.7, label='Theoretical $\\sqrt{N}$ improvement')

plt.xlabel('Bin Size')
plt.ylabel('Standard Deviation')
plt.title('Noise Reduction with Temporal Binning')
plt.legend()
plt.grid(True)
plt.show()

print("\nTemporal Binning Analysis:")
print(f"{'Bin Size':>10} {'Std Dev':>12}")
print("-" * 25)
for b, std in zip(bin_sizes, binned_stds):
    print(f"{b:10d} {std:12.6f}")

# differential light curve for mismatch case delta t = 2000 K

result_mismatch = differential_photometry(6000, 4000, "lsst2023-g",
                                          n_exposures, exposure_time, starting_photons)
time = np.arange(n_exposures)

plt.figure(figsize=(10, 6))

plt.plot(time, result_mismatch['diff_lc'], linewidth=0.5)
plt.axhline(y=np.mean(result_mismatch['diff_lc']), linestyle='--', linewidth=1,
            label=f'Mean = {np.mean(result_mismatch["diff_lc"]):.4f}')

plt.xlabel('Exposure Number')
plt.ylabel('Target / Reference')
plt.title('Differential Light Curve ($\Delta$T = 2000 K)')
plt.legend()
plt.grid(True)
plt.show()

print(f"Differential mean: {np.mean(result_mismatch['diff_lc']):.4f}")
print(f"Differential std: {np.std(result_mismatch['diff_lc']):.4f}")
print(f"Relative error: {np.std(result_mismatch['diff_lc'])/np.mean(result_mismatch['diff_lc']):.5f}")

# lsst filter curves

plt.figure(figsize=(10, 6))

for filt_name in filter_names:
    plt.plot(wavelengths_nm, filters_dict[filt_name], linewidth=2, label=filt_name[-1])

plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission')
plt.title('LSST Filter Response Curves')
plt.legend()
plt.grid(True)
plt.show()

# pipeline demo for 6000 K star in g band

T_demo = 6000

spectrum = planck(wavelengths_m, T_demo)
spectrum_norm = spectrum / np.trapezoid(spectrum, wavelengths_m)

transmission = atmosphere_transmission(wavelengths_nm, airmass=1.5)
after_atmosphere = spectrum_norm * transmission

after_filter = after_atmosphere * filters_dict["lsst2023-g"]

after_qe = after_filter * qe

plt.figure(figsize=(12, 6))

plt.plot(wavelengths_nm, spectrum_norm / np.max(spectrum_norm), linewidth=2, label='Blackbody')
plt.plot(wavelengths_nm, after_atmosphere / np.max(spectrum_norm), linewidth=2, label='After Atmosphere')
plt.plot(wavelengths_nm, after_filter / np.max(spectrum_norm), linewidth=2, label='After Filter (g-band)')
plt.plot(wavelengths_nm, after_qe / np.max(spectrum_norm), linewidth=2, label='After Detector QE')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Relative Flux')
plt.title('Simulation Pipeline: 6000K Star through g-band')
plt.legend()
plt.grid(True)
plt.xlim([300, 1000])
plt.show()

# qe curve

plt.figure(figsize=(10, 6))

plt.plot(wavelengths_nm, qe * 100, linewidth=2)

plt.xlabel('Wavelength (nm)')
plt.ylabel('Quantum Efficiency (%)')
plt.title('Detector Quantum Efficiency')
plt.grid(True)
plt.ylim([0, 100])
plt.show()

# cloud transparency variations

flux_test = run_pipeline(6000, "lsst2023-g", True)
expected_counts = flux_test * starting_photons * exposure_time
light_curve_data, cloud_factor = light_curve(expected_counts, n_exposures)
time = np.arange(n_exposures)

plt.figure(figsize=(10, 6))

plt.plot(time, cloud_factor, linewidth=0.5)
plt.axhline(y=0.90, linestyle='--', linewidth=2, label='Mean (0.90)')

plt.xlabel('Exposure Number')
plt.ylabel('Transparency Factor')
plt.title('Cloud Transparency Variations')
plt.legend()
plt.grid(True)
plt.ylim([0.7, 1.05])
plt.show()

# filter comparison bar chart

temps_plot = [4000, 6000, 10000]
fluxes_by_temp = {}

for T in temps_plot:
    fluxes = []
    for filt in filter_names:
        flux = run_pipeline(T, filt, True)
        fluxes.append(flux)
    fluxes_by_temp[T] = fluxes

x = np.arange(len(filter_names))
width = 0.25

plt.figure(figsize=(10, 6))

for i, T in enumerate(temps_plot):
    plt.bar(x + i*width, fluxes_by_temp[T], width, label=f'{T} K', alpha=0.8)

plt.xlabel('Filter')
plt.ylabel('Detected Flux (normalised)')
plt.title('Flux Through Different Filters')
plt.xticks(x + width, ['u', 'g', 'r', 'i', 'z'])
plt.legend()
plt.grid(True, alpha=0.3, axis='y')
plt.show()

# atmospheric loss vs temp

temps_loss = np.arange(3000, 12001, 1000)

plt.figure(figsize=(10, 6))

for filt in filter_names:
    losses = []
    for T in temps_loss:
        flux_no_atm = run_pipeline(T, filt, False)
        flux_atm = run_pipeline(T, filt, True)
        loss = 100 * (flux_no_atm - flux_atm) / flux_no_atm
        losses.append(loss)
    plt.plot(temps_loss, losses, 'o-', linewidth=2, markersize=6, label=filt[-1])

plt.xlabel('Temperature (K)')
plt.ylabel('Flux Loss (%)')
plt.title('Atmospheric Flux Loss vs Stellar Temperature')
plt.legend()
plt.grid(True)
plt.show()

# space vs ground comparison

filters_sg = ["lsst2023-g", "lsst2023-r", "lsst2023-i"]
delta_ts_sg = [0, 1000, 2000, 3000, 4000]

ground_errors = {f: [] for f in filters_sg}
space_errors = {f: [] for f in filters_sg}

for filt in filters_sg:
    for delta_t in delta_ts_sg:
        ref_temp = 6000 - delta_t
        if ref_temp > 0:
            # ground-based

            result_g = differential_photometry(6000, ref_temp, filt,
                                              n_exposures, exposure_time, starting_photons,
                                              airmass=1.5)
            ground_errors[filt].append(np.std(result_g['diff_lc']) / np.mean(result_g['diff_lc']))

            # space-based

            target_flux_s = run_pipeline(6000, filt, atmosphere=False)
            ref_flux_s = run_pipeline(ref_temp, filt, atmosphere=False)

            target_counts_s = target_flux_s * starting_photons * exposure_time
            ref_counts_s = ref_flux_s * starting_photons * exposure_time

            target_lc_s = add_noise(np.full(n_exposures, target_counts_s))
            ref_lc_s = add_noise(np.full(n_exposures, ref_counts_s))
            diff_lc_s = target_lc_s / ref_lc_s

            space_errors[filt].append(np.std(diff_lc_s) / np.mean(diff_lc_s))

plt.figure(figsize=(10, 6))

for filt in filters_sg:
    plt.plot(delta_ts_sg, ground_errors[filt], 'o-', linewidth=2, markersize=6, label=f'{filt[-1]} Ground')
    plt.plot(delta_ts_sg, space_errors[filt], 's--', linewidth=2, markersize=6, label=f'{filt[-1]} Space')

plt.xlabel('Temperature Difference (K)')
plt.ylabel('Relative Error')
plt.title('Ground vs Space Differential Photometry')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

print("\nGround vs Space Comparison:")
print(f"{'Filter':>8} {'ΔT (K)':>8} {'Ground Error':>14} {'Space Error':>14}")
print("-" * 48)
for filt in filters_sg:
    for i, dt in enumerate(delta_ts_sg):
        print(f"{filt[-1]:>8} {dt:8d} {ground_errors[filt][i]:14.5f} {space_errors[filt][i]:14.5f}")

# wien's law validation

temperatures_wien = [3000, 4000, 5000, 6000, 8000, 10000]
peak_wavelengths = [wien_peak(T) for T in temperatures_wien]
theory_peaks = [2.898e6 / T for T in temperatures_wien]

plt.figure(figsize=(8, 6))

plt.plot(temperatures_wien, peak_wavelengths, 'o-', linewidth=2, markersize=8, label='Calculated')
plt.plot(temperatures_wien, theory_peaks, '--', linewidth=2, alpha=0.7, label="Wien's Law (theory)")

plt.xlabel('Temperature (K)')
plt.ylabel('Peak Wavelength (nm)')
plt.title("Wien's Displacement Law Validation")
plt.legend()
plt.grid(True)
plt.show()

print("Wien's Law Validation:")
print(f"{'Temp (K)':>10} {'Peak (nm)':>10} {'Theory (nm)':>12}")
for T, peak, theory in zip(temperatures_wien, peak_wavelengths, theory_peaks):
    print(f"{T:10d} {peak:10.1f} {theory:12.1f}")

# summary of results

print("FULL RESULTS SUMMARY")

print("\nSpectral Mismatch Analysis (Target: 6000K, Filter: lsst2023-g):")
print(f"{'ΔT':>8} {'Mean Ratio':>12} {'Std Dev':>12} {'Rel Error':>12}")

delta_t_range = np.arange(-4000, 4001, 500)
for delta_t in delta_t_range:
    ref_temp = 6000 - delta_t
    if ref_temp > 0:
        result = differential_photometry(6000, ref_temp, "lsst2023-g",
                                        n_exposures, exposure_time, starting_photons)
        mean_r = np.mean(result['diff_lc'])
        std_r = np.std(result['diff_lc'])
        rel_e = std_r / mean_r
        print(f"{delta_t:8d} {mean_r:12.4f} {std_r:12.4f} {rel_e:12.5f}")

print("\nTemporal Binning Analysis:")
print(f"{'Bin Size':>10} {'Std Dev':>12} {'Noise Reduction':>18}")

result_bin = differential_photometry(6000, 4000, "lsst2023-g",
                                    n_exposures, exposure_time, starting_photons)

bin_sizes = [1, 2, 5, 10, 20, 50]
base_std = np.std(result_bin['diff_lc'])
for b in bin_sizes:
    binned = bin_data(result_bin['diff_lc'], b)
    binned_std = np.std(binned)
    reduction = (1 - binned_std/base_std) * 100
    print(f"{b:10d} {binned_std:12.6f} {reduction:17.1f}%")

print("\nFilter Comparison:")
print(f"{'Filter':>8} {'ΔT=1000K':>12} {'ΔT=2000K':>12} {'ΔT=3000K':>12} {'ΔT=4000K':>12}")

for filt in filter_names:
    errors = []
    for delta_t in [1000, 2000, 3000, 4000]:
        ref_temp = 6000 - delta_t
        result = differential_photometry(6000, ref_temp, filt,
                                        n_exposures, exposure_time, starting_photons)
        errors.append(np.std(result['diff_lc']) / np.mean(result['diff_lc']))
    print(f"{filt[-1]:>8} {errors[0]:12.5f} {errors[1]:12.5f} {errors[2]:12.5f} {errors[3]:12.5f}")

# blackbody spectra with filter bandpass

plt.figure(figsize=(12, 6))

for T in [4000, 6000, 10000]:
    spectrum = planck(wavelengths_m, T)
    spectrum_norm = spectrum / np.max(spectrum)
    plt.plot(wavelengths_nm, spectrum_norm, linewidth=2, label=f'{T} K')

plt.fill_between(wavelengths_nm, 0, filters_dict["lsst2023-u"] * 0.8, alpha=0.3, label='u-band')
plt.fill_between(wavelengths_nm, 0, filters_dict["lsst2023-i"] * 0.8, alpha=0.3, label='i-band')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalised Flux / Filter Response')
plt.title('Stellar Spectra with Filter Bandpasses')
plt.legend()
plt.grid(True)
plt.xlim([300, 1100])
plt.show()

# distribution comparison delta t = 0 and 2000

result_match = differential_photometry(6000, 6000, "lsst2023-g",
                                       n_exposures, exposure_time, starting_photons)
result_mismatch = differential_photometry(6000, 4000, "lsst2023-g",
                                          n_exposures, exposure_time, starting_photons)

plt.figure(figsize=(10, 6))

plt.hist(result_match['diff_lc'], bins=30, alpha=0.5, label='$\Delta$T = 0 K')
plt.hist(result_mismatch['diff_lc'], bins=30, alpha=0.5, label='$\Delta$T = 2000 K')
plt.axvline(x=1.0, linestyle='--', linewidth=2, label='Ratio = 1.0')

plt.xlabel('Target / Reference Ratio')
plt.ylabel('Frequency')
plt.title('Distribution Comparison: Matched vs Mismatched Temperatures')
plt.legend()
plt.grid(True)
plt.show()

# multiple reference stars

def multiple_reference_differential_photometry(target_temp, reference_temps, filter_name,
                                               n_exposures=350, exposure_time=10,
                                               starting_photons=100000, airmass=1.5):

    target_flux = run_pipeline(target_temp, filter_name, True, airmass)
    target_counts = target_flux * starting_photons * exposure_time

    cloud_factor = np.random.normal(0.90, 0.03, n_exposures)
    cloud_factor = np.clip(cloud_factor, 0, 1)

    all_ref_lcs = []
    ref_temps_used = []

    for ref_temp in reference_temps:
        ref_flux = run_pipeline(ref_temp, filter_name, True, airmass)
        ref_counts = ref_flux * starting_photons * exposure_time

        ref_lc = []
        for transparency in cloud_factor:
            r_measured = add_noise(ref_counts * transparency)
            ref_lc.append(r_measured)

        all_ref_lcs.append(np.array(ref_lc))
        ref_temps_used.append(ref_temp)

    avg_ref_lc = np.mean(all_ref_lcs, axis=0)

    target_lc = []
    for transparency in cloud_factor:
        t_measured = add_noise(target_counts * transparency)
        target_lc.append(t_measured)

    target_lc = np.array(target_lc)

    diff_lc = target_lc / avg_ref_lc

    return {
        'target_lc': target_lc,
        'avg_ref_lc': avg_ref_lc,
        'diff_lc': diff_lc,
        'ref_temps': ref_temps_used,
        'target_temp': target_temp
    }

target_temp = 6000

result_single = differential_photometry(target_temp, 4000, "lsst2023-g",
                                        n_exposures, exposure_time, starting_photons)
single_error = np.std(result_single['diff_lc']) / np.mean(result_single['diff_lc'])

result_two = multiple_reference_differential_photometry(
    target_temp, [5500, 4000], "lsst2023-g",
    n_exposures, exposure_time, starting_photons)
two_error = np.std(result_two['diff_lc']) / np.mean(result_two['diff_lc'])

result_three = multiple_reference_differential_photometry(
    target_temp, [6500, 5500, 4000], "lsst2023-g",
    n_exposures, exposure_time, starting_photons)
three_error = np.std(result_three['diff_lc']) / np.mean(result_three['diff_lc'])

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(np.arange(n_exposures), result_single['target_lc'], linewidth=0.5, label='Target (6000K)')
axes[0].plot(np.arange(n_exposures), result_single['reference_lc'], linewidth=0.5, label='Ref (4000K)')
axes[0].set_xlabel('Exposure Number')
axes[0].set_ylabel('Detected Counts')
axes[0].set_title('Single Reference Star')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

configs = ['Single\n(4000K)', 'Two\n(5500+4000K)', 'Three\n(6500+5500+4000K)']
errors = [single_error, two_error, three_error]

axes[1].bar(configs, errors, alpha=0.7, edgecolor='black')
for i, (config, err) in enumerate(zip(configs, errors)):
    axes[1].text(i, err, f'{err:.5f}', ha='center', va='bottom')
axes[1].set_ylabel('Relative Error')
axes[1].set_title('Multiple Reference Stars Comparison')
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

print("Multiple Reference Stars Results:")
print(f"  Single reference (4000K):           {single_error:.5f}")
print(f"  Two references (5500K, 4000K):      {two_error:.5f}")
print(f"  Three references (6500, 5500, 4000): {three_error:.5f}")

# flux change analysis

print(f"\n{'ΔT':>6} {'Target Flux':>14} {'Ref Flux':>14} {'Flux Ratio':>12} {'Count Ratio':>14}")
print("-" * 65)

for delta_t in range(0, 4001, 1000):
    ref_temp = 6000 - delta_t
    if ref_temp <= 0:
        continue

    target_flux = run_pipeline(6000, "lsst2023-g", True)
    ref_flux = run_pipeline(ref_temp, "lsst2023-g", True)

    target_counts = target_flux * starting_photons * exposure_time
    ref_counts = ref_flux * starting_photons * exposure_time

    flux_ratio = target_flux / ref_flux
    count_ratio = target_counts / ref_counts

    print(f"{delta_t:6d} {target_flux:14.6f} {ref_flux:14.6f} {flux_ratio:12.4f} {count_ratio:14.4f}")

print(f"\n{'Temp (K)':>10} {'No Atm Flux':>14} {'With Atm Flux':>14} {'Loss (%)':>10}")
print("-" * 52)

for T in [4000, 6000, 10000]:
    flux_no_atm = run_pipeline(T, "lsst2023-g", False)
    flux_with_atm = run_pipeline(T, "lsst2023-g", True)
    loss = 100 * (flux_no_atm - flux_with_atm) / flux_no_atm
    print(f"{T:10d} {flux_no_atm:14.6f} {flux_with_atm:14.6f} {loss:10.2f}")

# snr analysis

print(f"\n{'ΔT':>6} {'Target SNR':>12} {'Ref SNR':>12} {'Diff SNR':>12} {'SNR Ratio':>12}")
print("-" * 60)

for delta_t in range(0, 4001, 1000):
    ref_temp = 6000 - delta_t
    if ref_temp <= 0:
        continue

    result = differential_photometry(6000, ref_temp, "lsst2023-g",
                                    n_exposures, exposure_time, starting_photons)

    target_snr = np.mean(result['target_lc']) / np.std(result['target_lc'])
    ref_snr = np.mean(result['reference_lc']) / np.std(result['reference_lc'])
    diff_snr = np.mean(result['diff_lc']) / np.std(result['diff_lc'])
    snr_ratio = diff_snr / target_snr

    print(f"{delta_t:6d} {target_snr:12.2f} {ref_snr:12.2f} {diff_snr:12.2f} {snr_ratio:12.2f}")

print(f"\n{'Exp (s)':>8} {'Target SNR':>12} {'Diff SNR':>12}")
print("-" * 36)

for et in [1, 5, 10, 30, 60]:
    result = differential_photometry(6000, 4000, "lsst2023-g",
                                    n_exposures, exposure_time=et,
                                    starting_photons=starting_photons)

    target_snr = np.mean(result['target_lc']) / np.std(result['target_lc'])
    diff_snr = np.mean(result['diff_lc']) / np.std(result['diff_lc'])

    print(f"{et:8d} {target_snr:12.2f} {diff_snr:12.2f}")

# source brightness effect

brightness_levels = [10000, 50000, 100000, 500000, 1000000]
brightness_errors = []
brightness_snrs = []

for photons in brightness_levels:
    result = differential_photometry(6000, 4000, "lsst2023-g",
                                    n_exposures, exposure_time,
                                    starting_photons=photons)

    rel_error = np.std(result['diff_lc']) / np.mean(result['diff_lc'])
    diff_snr = np.mean(result['diff_lc']) / np.std(result['diff_lc'])

    brightness_errors.append(rel_error)
    brightness_snrs.append(diff_snr)

plt.figure(figsize=(10, 6))

plt.plot(brightness_levels, brightness_errors, 'o-', linewidth=2, markersize=8)

plt.xlabel('Source Brightness (photons)')
plt.ylabel('Relative Error')
plt.title('Effect of Source Brightness on Spectral Mismatch ($\Delta$T = 2000 K, g-band)')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.show()
