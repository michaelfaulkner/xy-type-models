from scipy import signal
import importlib
import numpy as np
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_modules = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_modules)
sample_getter = importlib.import_module("sample_getter")


def get_power_spectrum(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                       no_of_equilibration_sweeps, sampling_frequency=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    mean_zero_time_series = get_mean_zero_time_series(power_spectrum_string, output_directory, temperature_directory,
                                                      beta, no_of_sites, no_of_equilibration_sweeps)
    return get_normalised_component_averaged_power_spectrum(mean_zero_time_series, sampling_frequency)


def get_correlator_power_spectrum(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps, sampling_frequency=None, shifted_time_period=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    mean_zero_time_series = get_mean_zero_time_series(power_spectrum_string, output_directory, temperature_directory,
                                                      beta, no_of_sites, no_of_equilibration_sweeps)
    if shifted_time_period is None:
        shifted_time_period = int(0.1 * len(mean_zero_time_series[0]))
    shifted_mean_zero_time_series = np.roll(mean_zero_time_series, shifted_time_period, axis=1)
    two_point_correlator = np.conj(mean_zero_time_series) * shifted_mean_zero_time_series
    return get_normalised_component_averaged_power_spectrum(two_point_correlator, sampling_frequency)


def get_sampling_frequency(output_directory, sampling_frequency, temperature_directory):
    if sampling_frequency is None:
        acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        physical_time_scale = acceptance_rates[1] * acceptance_rates[0] ** 2 / 24.0
        sampling_frequency = 1.0 / physical_time_scale
    return sampling_frequency


def get_mean_zero_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                              no_of_equilibration_sweeps):
    get_sample_method = getattr(sample_getter, "get_" + power_spectrum_string)
    sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)[no_of_equilibration_sweeps:]
    mean_zero_sample = sample - np.mean(sample, axis=0)
    mean_zero_sample = np.atleast_2d(mean_zero_sample)
    if len(mean_zero_sample) > 1:
        mean_zero_sample = mean_zero_sample.transpose()
    return mean_zero_sample


def get_normalised_component_averaged_power_spectrum(mean_zero_time_series, sampling_frequency):
    power_spectra = np.atleast_2d(
        [signal.periodogram(component, fs=sampling_frequency) for component in mean_zero_time_series])
    # now average over Cartesian components of original sample, where the `[:, 1:]' removes the f = 0 value as...
    # ...this value is invalid for a finite-time signal
    component_averaged_power_spectrum = np.mean(power_spectra, axis=0)[:, 1:]
    # normalise power spectrum with respect to its low-frequency value
    return (component_averaged_power_spectrum[0],
            component_averaged_power_spectrum[1] / component_averaged_power_spectrum[1, 0])


'''
def get_magnetisation_phase_correlator_power_spectrum(output_directory, temperature_directory, beta, no_of_sites,
                                                      no_of_equilibration_sweeps, sampling_frequency=None,
                                                      shifted_time_period=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    sample = sample_getter.get_magnetisation_phase(output_directory, temperature_directory, beta, no_of_sites)[
             no_of_equilibration_sweeps:]
    relative_sample = sample - np.mean(sample)
    if shifted_time_period is None:
        shifted_time_period = int(0.1 * len(sample))
    shifted_relative_sample = np.roll(relative_sample, shifted_time_period)
    two_point_correlator = np.conj(relative_sample) * shifted_relative_sample
    power_spectrum = signal.periodogram(two_point_correlator, fs=sampling_frequency)
    return power_spectrum[0], np.array([item if abs(item) > 1.0e-15 else 0.0 for item in power_spectrum[1]])
'''
