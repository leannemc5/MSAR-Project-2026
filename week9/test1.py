# compare spectral mismatch effects across different photometric filters

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k
from scipy.interpolate import interp1d
from speclite import filters

# input parameters
starting_photons = 100000
exposure_time = 10
n_exposures = 350
read_noise_sigma = 100

target_temperature = 6000

wavelength_min = 300
wavelength_max = 1100
n_wavelengths = 5000

filter_names = ["lsst2023-u", "lsst2023-g", "lsst2023-r", "lsst2023-i", "lsst2023-z"]

# qe data
qe_wavelength_nm = np.array([350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000])

qe_percent = np.array([20, 39, 55, 63, 72, 77, 85, 87, 86, 85, 80, 60, 38, 17])

# wavelength grid
wavelengths_nm = np.linspace(wavelength_min, wavelength_max, n_wavelengths)
wavelengths_m = wavelengths_nm * 1e-9

# planck function
def planck(wavelength_m, temperature):

    exponent = (h * c) / (wavelength_m * k * temperature)

    exponent = np.clip(exponent, 0, 700)

    intensity = (2 * h * c**2) / (wavelength_m**5 * (np.exp(exponent) - 1))

    return intensity

# atmospheric extinction
def atmosphere_transmission(wavelengths_nm, airmass=1.0):

    tau_550 = 0.1

    tau_rayleigh = tau_550 * (550 / wavelengths_nm)**4

    tau_total = tau_rayleigh * airmass

    transmission = np.exp(-tau_total)

    return transmission

# detector qe
def qe_curve(wavelengths_nm):

    qe_interp = interp1d(qe_wavelength_nm, qe_percent / 100, bounds_error=False, fill_value=0)

    return qe_interp(wavelengths_nm)

# load filters
def filter_responses(wavelengths_nm):

    responses = {}

    for filt_name in filter_names:

        filt = filters.load_filter(filt_name)

        response_func = interp1d(filt.wavelength / 10, filt.response, bounds_error=False, fill_value=0)

        responses[filt_name] = response_func(wavelengths_nm)

    return responses

qe = qe_curve(wavelengths_nm)
filters_dict = filter_responses(wavelengths_nm)

# pipeline function
def run_pipeline(temperature, filter_name, atmosphere=True, airmass=1.5):

    spectrum = planck(wavelengths_m, temperature)

    spectrum /= np.trapezoid(spectrum, wavelengths_m)

    if atmosphere:

        transmission = atmosphere_transmission(wavelengths_nm, airmass)

        spectrum *= transmission

    spectrum *= filters_dict[filter_name]

    spectrum *= qe

    detected_flux = np.trapezoid(spectrum, wavelengths_m)

    return detected_flux

# noise model
def add_noise(counts, read_noise_sigma=100):

    poisson_counts = np.random.poisson(np.maximum(counts, 0))

    read_noise = np.random.normal(0, read_noise_sigma, size=np.shape(poisson_counts))

    return poisson_counts + read_noise

# differential photometry
def differential_photometry(target_temp, reference_temp, filter_name, n_exposures=350, exposure_time=10, starting_photons=100000, airmass=1.5):

    target_flux = run_pipeline(target_temp, filter_name, True, airmass)

    reference_flux = run_pipeline(reference_temp, filter_name, True, airmass)

    target_counts = target_flux * starting_photons * exposure_time

    reference_counts = reference_flux * starting_photons * exposure_time

    cloud_factor = np.random.normal(0.90, 0.03, n_exposures)

    cloud_factor = np.clip(cloud_factor, 0, 1)

    target_lc = []
    reference_lc = []

    for transparency in cloud_factor:

        t_measured = add_noise(target_counts * transparency)
        r_measured = add_noise(reference_counts * transparency)

        target_lc.append(t_measured)
        reference_lc.append(r_measured)

    target_lc = np.array(target_lc)
    reference_lc = np.array(reference_lc)

    diff_lc = target_lc / reference_lc

    return {
        'target_lc': target_lc,
        'reference_lc': reference_lc,
        'diff_lc': diff_lc,
        'target_flux': target_flux,
        'reference_flux': reference_flux,
        'target_counts': target_counts,
        'reference_counts': reference_counts
    }

# filter comparison
delta_ts_test = [1000, 2000, 3000, 4000]

filter_errors = {}

for filt in filter_names:

    errors = []

    for delta_t in delta_ts_test:

        ref_temp = target_temperature - delta_t

        if ref_temp > 0:

            result = differential_photometry(target_temperature, ref_temp, filt, n_exposures, exposure_time, starting_photons)

            rel_error = np.std(result['diff_lc']) / np.mean(result['diff_lc'])

            errors.append(rel_error)

    filter_errors[filt] = errors

print()
print("filter comparison")
print(f"{'Filter':>8} {'ΔT=1000K':>12} {'ΔT=2000K':>12} {'ΔT=3000K':>12} {'ΔT=4000K':>12}")
print("-" * 56)

for filt in filter_names:

    print(f"{filt[-1]:>8} {filter_errors[filt][0]:12.5f} {filter_errors[filt][1]:12.5f} {filter_errors[filt][2]:12.5f} {filter_errors[filt][3]:12.5f}")

# plot
plt.figure(figsize=(10,6))

for filt in filter_names:

    plt.plot(delta_ts_test, filter_errors[filt], 'o-', linewidth=2, markersize=6, label=filt[-1])

plt.xlabel("Temperature Difference (K)")
plt.ylabel("Relative Error")
plt.title("Spectral Mismatch Across Different Filters")
plt.legend()
plt.grid(True)
plt.show()
