from scipy import signal
import importlib
import math
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
                                     no_of_equilibration_sweeps, time_period_shift=10, sampling_frequency=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    if time_period_shift >= len(time_series[0]):
        raise Exception("time_period_shift must be an integer less than the sample size.")
    return get_component_averaged_power_spectrum(
        get_two_point_correlator(time_series - np.mean(time_series, axis=1), time_period_shift), sampling_frequency)


def get_power_trispectrum(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                          no_of_equilibration_sweeps, no_of_jobs, pool, no_of_octaves=2, base_time_period_shift=1,
                          sampling_frequency=None):
    if no_of_octaves <= 0:
        raise Exception("no_of_octaves must be a positive integer.")
    if no_of_jobs == 1:
        power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(power_spectrum_string,
                                                                                    output_directory,
                                                                                    temperature_directory, beta,
                                                                                    no_of_sites,
                                                                                    no_of_equilibration_sweeps,
                                                                                    base_time_period_shift,
                                                                                    no_of_octaves, sampling_frequency)
    else:
        power_spectra_of_correlators = pool.starmap(get_power_spectra_of_trispectrum_correlators,
                                                    [(power_spectrum_string, f"{output_directory}/job_{job_number + 1}",
                                                      temperature_directory, beta, no_of_sites,
                                                      no_of_equilibration_sweeps, base_time_period_shift, no_of_octaves,
                                                      sampling_frequency)
                                                     for job_number in range(no_of_jobs)])
        power_spectra_of_correlators = np.mean(np.array(power_spectra_of_correlators, dtype=object), axis=0)
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    spectra_in_frequency_shift_space = [np.absolute(item) for index, item in
                                        enumerate(np.fft.fft(transposed_power_spectra).transpose())
                                        if (index == 0 or index == 2 ** (math.floor(math.log(index, 2))))]
    # normalise power trispectrum with respect to its low-frequency value
    spectra_in_frequency_shift_space = [spectrum / spectrum[0] for spectrum in spectra_in_frequency_shift_space]
    return [np.atleast_1d(np.fft.fftfreq(len(transposed_power_spectra[0]), d=base_time_period_shift)[1]),
            power_spectra_of_correlators[0, 0], spectra_in_frequency_shift_space]


def get_power_trispectrum_estimator(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                    no_of_equilibration_sweeps, no_of_jobs, pool, no_of_octaves=2,
                                    base_time_period_shift=1, sampling_frequency=None):
    if no_of_octaves <= 0:
        raise Exception("no_of_octaves must be a positive integer.")
    if no_of_jobs == 1:
        power_trispectrum_estimator = get_single_observation_of_power_trispectrum_estimator(power_spectrum_string,
                                                                                            output_directory,
                                                                                            temperature_directory, beta,
                                                                                            no_of_sites,
                                                                                            no_of_equilibration_sweeps,
                                                                                            no_of_octaves,
                                                                                            base_time_period_shift,
                                                                                            sampling_frequency)
    else:
        power_trispectrum_estimators = pool.starmap(get_single_observation_of_power_trispectrum_estimator,
                                                    [(power_spectrum_string, f"{output_directory}/job_{job_number + 1}",
                                                      temperature_directory, beta, no_of_sites,
                                                      no_of_equilibration_sweeps, no_of_octaves, base_time_period_shift,
                                                      sampling_frequency)
                                                     for job_number in range(no_of_jobs)])
        power_trispectrum_estimator = np.mean(np.array(power_trispectrum_estimators, dtype=object), axis=0)
    # normalise estimator of power trispectrum with respect to its low-frequency value
    power_trispectrum_estimator[2] = [spectrum / spectrum[0] for spectrum in power_trispectrum_estimator[2]]
    return power_trispectrum_estimator


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
        time series (rather than using np.roll()).  Previously, we returned
        np.conj(time_series[:, time_period_shift:]) * time_series[:, :len(time_series[0]) - time_period_shift],
        but have since removed the np.conj() operation as we only consider real-valued signals."""
    return time_series[:, time_period_shift:] * time_series[:, :len(time_series[0]) - time_period_shift]


def get_power_spectra_of_trispectrum_correlators(power_spectrum_string, output_directory, temperature_directory, beta,
                                                 no_of_sites, no_of_equilibration_sweeps, base_time_period_shift,
                                                 no_of_octaves, sampling_frequency):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(power_spectrum_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    if base_time_period_shift * 2 ** (no_of_octaves + 1) >= len(time_series[0]):
        raise Exception("base_time_period_shift * 2 ** (no_of_octaves + 1) must be less than the smallest power of 2 "
                        "that is less than the sample size.")
    # create 2 ** (no_of_octaves + 1) two-point correlators; here, [:, :len(time_series[0]) -
    # (2 ** (no_of_octaves + 1) - 1) * base_time_period_shift] ensures that (each component of) all correlators are the
    # same length, which ensures that their power spectra have common frequency values
    correlators = [
        get_two_point_correlator(time_series - np.mean(time_series, axis=1), i * base_time_period_shift)[
            :, :len(time_series[0]) - (2 ** (no_of_octaves + 1) - 1) * base_time_period_shift]
        for i in range(2 ** (no_of_octaves + 1))]
    return np.array(
        [get_component_averaged_power_spectrum(correlator, sampling_frequency) for correlator in correlators])


def get_single_observation_of_power_trispectrum_estimator(power_spectrum_string, output_directory,
                                                          temperature_directory, beta, no_of_sites,
                                                          no_of_equilibration_sweeps, no_of_octaves=2,
                                                          base_time_period_shift=1, sampling_frequency=None):
    power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(power_spectrum_string, output_directory,
                                                                                temperature_directory, beta,
                                                                                no_of_sites, no_of_equilibration_sweeps,
                                                                                base_time_period_shift, no_of_octaves,
                                                                                sampling_frequency)
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    norm_of_spectra_in_frequency_shift_space = np.array(
        [np.absolute(item) for index, item in enumerate(np.fft.fft(transposed_power_spectra).transpose())
         if (index == 0 or index == 2 ** (math.floor(math.log(index, 2))))])
    return [np.atleast_1d(np.fft.fftfreq(len(transposed_power_spectra[0]), d=base_time_period_shift)[1]),
            power_spectra_of_correlators[0, 0], norm_of_spectra_in_frequency_shift_space]


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
