import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k
from scipy.interpolate import interp1d
from speclite import filters

# constants
starting_photons = 100000
exposure_time = 10
n_exposures = 350
read_noise_sigma = 100

wavelength_min = 300
wavelength_max = 1100
n_wavelengths = 5000

filter_names = ["lsst2023-u", "lsst2023-g", "lsst2023-r", "lsst2023-i", "lsst2023-z"]

# wavelength grid
wavelengths_nm = np.linspace(wavelength_min, wavelength_max, n_wavelengths)
wavelengths_m = wavelengths_nm * 1e-9

# planck function
def planck(wavelength_m, temperature):

    exponent = (h * c) / (wavelength_m * k * temperature)

    exponent = np.clip(exponent, 0, 700)

    intensity = (2 * h * c**2) / (wavelength_m**5 * (np.exp(exponent) - 1))

    return intensity

# wien's law validation
def wien_peak(temperature):

    b = 2.898e-3

    return (b / temperature) * 1e9

# atmospheric extinction
def atmosphere_transmission(wavelengths_nm, airmass=1.0):

    tau_550 = 0.1

    tau_rayleigh = tau_550 * (550 / wavelengths_nm)**4

    tau_total = tau_rayleigh * airmass

    transmission = np.exp(-tau_total)

    return transmission

# detector quantum efficiency
def qe_curve(wavelengths_nm):

    qe_wavelength_nm = np.array([350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000])

    qe_percent = np.array([20, 39, 55, 63, 72, 77, 85, 87, 86, 85, 80, 60, 38, 17])

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

# pipeline
def run_pipeline(temperature, filter_name, atmosphere=True, airmass=1.0):

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

# lightcurve with clouds
def light_curve(expected_counts, n_exposures, cloud_std=0.03, read_noise_sigma=100):

    cloud_factor = np.random.normal(0.90, cloud_std, n_exposures)

    cloud_factor = np.clip(cloud_factor, 0, 1)

    curve = []

    for transparency in cloud_factor:

        counts = expected_counts * transparency

        noisy_counts = add_noise(counts, read_noise_sigma)

        curve.append(noisy_counts)

    return np.array(curve), cloud_factor

# temporal binning
def bin_data(data, bin_size):

    n = (len(data) // bin_size) * bin_size

    return data[:n].reshape(-1, bin_size).mean(axis=1)

# differential photometry
def differential_photometry(target_temp, reference_temp, filter_name, n_exposures=350, exposure_time=10, starting_photons=100000, airmass=1.0):

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
        'reference_counts': reference_counts}
