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
    get_sample_method = getattr(sample_getter, "get_" + power_spectrum_string)
    sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)[no_of_equilibration_sweeps:]
    # sample_variance = np.sum(np.var(np.atleast_2d(sample), axis=1))
    # the following line subtracts the Gaussian contribution ot the power spectrum,
    # though numerical experiments seem to indicate it's redundant
    sample -= np.mean(sample, axis=0)
    sample = np.atleast_2d(sample)
    if len(sample) > 1:
        sample = sample.transpose()
    power_spectrum = np.atleast_2d([signal.periodogram(component, fs=sampling_frequency) for component in sample])
    component_averaged_power_spectrum = np.mean(power_spectrum, axis=0)
    return component_averaged_power_spectrum[0], np.array([item if abs(item) > 1.0e-15 else 0.0 for item in
                                                           component_averaged_power_spectrum[1]])


def get_correlator_power_spectrum(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps, sampling_frequency=None, shifted_time_period=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    get_sample_method = getattr(sample_getter, "get_" + power_spectrum_string)
    sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)[no_of_equilibration_sweeps:]
    relative_sample = sample - np.mean(sample)
    if shifted_time_period is None:
        shifted_time_period = int(0.1 * len(sample))
    shifted_relative_sample = np.roll(relative_sample, shifted_time_period)
    two_point_correlator = np.conj(relative_sample) * shifted_relative_sample
    return signal.periodogram(two_point_correlator, fs=sampling_frequency)


def get_sampling_frequency(output_directory, sampling_frequency, temperature_directory):
    if sampling_frequency is None:
        acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        physical_time_scale = acceptance_rates[1] * acceptance_rates[0] ** 2 / 24.0
        sampling_frequency = 1.0 / physical_time_scale
    return sampling_frequency