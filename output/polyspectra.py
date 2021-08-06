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
    time_series = get_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    return get_component_averaged_power_spectrum(time_series, sampling_frequency)


def get_power_spectrum_of_correlator(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                     no_of_equilibration_sweeps, sampling_frequency=None, time_period_shift=10):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    if time_period_shift >= len(time_series[0]):
        raise Exception("time_period_shift must be an integer less than the sample size.")
    return get_component_averaged_power_spectrum(
        get_two_point_correlator(time_series - np.mean(time_series, axis=0), time_period_shift), sampling_frequency)


def get_power_trispectrum(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                          no_of_equilibration_sweeps, sampling_frequency=None, no_of_octaves=2):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    if no_of_octaves <= 0:
        raise Exception("no_of_octaves must be a positive integer.")
    if 2 ** no_of_octaves >= len(time_series[0]):
        raise Exception("2 ** no_of_octaves must be less than the sample size.")
    # remove supplementary elements of each component i of time_series (w.r.t. enforcing
    # len(time_series[i]) % (2 ** no_of_octaves) == 0)
    if len(time_series[0]) % (2 ** no_of_octaves) != 0:
        time_series = time_series[:, :len(time_series[0]) - len(time_series[0]) // (2 ** no_of_octaves) - 1]
    base_time_period_shift = int(len(time_series[0]) / (2 ** no_of_octaves))
    # create 2 ** no_of_octaves two-point correlators
    correlators = [get_two_point_correlator(time_series - np.mean(time_series, axis=0), i * base_time_period_shift)
                   for i in range(2 ** no_of_octaves)]
    # in the following, [:].reshape(-1, 2 ** (no_of_octaves - index - 1)))[:, 0].reshape(2, int(len(correlator[0])
    # / (2 ** (no_of_octaves - index - 1))))) removes each element from each power_spectrum whose frequency value does
    # not correspond to one of the frequency values of the power spectrum of the shortest correlator
    power_spectra_of_correlators = np.array([
        get_component_averaged_power_spectrum(correlator, sampling_frequency)[:].reshape(-1, 2 ** (
                no_of_octaves - index - 1))[:, 0].reshape(2, int(len(correlator[0]) / (2 ** (
                    no_of_octaves - index - 1)))) for index, correlator in enumerate(correlators)])
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    return [np.fft.fftfreq(len(transposed_power_spectra[0]), d=base_time_period_shift),
            power_spectra_of_correlators[0, 0], np.fft.fft(transposed_power_spectra).transpose()]


def get_sampling_frequency(output_directory, sampling_frequency, temperature_directory):
    if sampling_frequency is None:
        acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        physical_time_scale = acceptance_rates[1] * acceptance_rates[0] ** 2 / 24.0
        sampling_frequency = 1.0 / physical_time_scale
    return sampling_frequency


def get_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                    no_of_equilibration_sweeps):
    get_sample_method = getattr(sample_getter, "get_" + power_spectrum_string)
    sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)[no_of_equilibration_sweeps:]
    sample = np.atleast_2d(sample)
    if len(sample) > 1:
        sample = sample.transpose()
    return sample


def get_component_averaged_power_spectrum(time_series, sampling_frequency):
    power_spectra = np.atleast_2d([signal.periodogram(component, fs=sampling_frequency) for component in
                                   time_series - np.mean(time_series, axis=1)])
    # now average over Cartesian components of original sample, where the `[:, 1:]' removes the f = 0 value as...
    # ...this value is invalid for a finite-time signal
    return np.mean(power_spectra, axis=0)[:, 1:]


def get_two_point_correlator(time_series, time_period_shift):
    """As the time series is not periodic, we chop off the first / last time_period_shift elements of each copy of the
        time series (rather than using np.roll())"""
    return np.conj(time_series[:, time_period_shift:]) * time_series[:, :len(time_series[0]) - time_period_shift]


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
