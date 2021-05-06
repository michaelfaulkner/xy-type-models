from scipy import signal
import importlib
import numpy as np
import os
import sys


# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.abspath(this_directory + '/../')
sys.path.insert(0, output_directory)
sample_getter = importlib.import_module('sample_getter')


def get_power_spectrum(sample, output_directory, temperature_directory, sampling_frequency=None):
    if sampling_frequency is None:
        acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        physical_time_scale = acceptance_rates[1] * acceptance_rates[0] ** 2 / 24.0
        sampling_frequency = 1.0 / physical_time_scale
    sample = np.atleast_2d(sample)
    if len(sample) > 1:
        sample = sample.transpose()
    power_spectrum = np.atleast_2d([signal.periodogram(component, fs=sampling_frequency) for component in sample])
    component_averaged_power_spectrum = np.mean(power_spectrum, axis=0)
    return component_averaged_power_spectrum[0], np.array([item if abs(item) > 1.0e-15 else 0.0 for item in
                                                           component_averaged_power_spectrum[1]])
