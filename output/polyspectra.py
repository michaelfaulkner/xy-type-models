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


def get_power_spectrum(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                       no_of_equilibration_sweeps, sampling_frequency=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    return get_component_averaged_power_spectrum(time_series, sampling_frequency)


def get_power_spectrum_of_correlator(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                                     no_of_equilibration_sweeps, time_period_shift=10, sampling_frequency=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                                  no_of_equilibration_sweeps)
    if time_period_shift >= len(time_series[0]):
        raise Exception("time_period_shift must be an integer less than the sample size.")
    return get_component_averaged_power_spectrum(
        get_two_point_correlator(time_series - np.mean(time_series, axis=1), time_period_shift), sampling_frequency)


def get_power_trispectrum(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                          no_of_equilibration_sweeps, no_of_jobs, pool, no_of_octaves=2, base_time_period_shift=1,
                          sampling_frequency=None):
    if no_of_octaves <= 0:
        raise Exception("no_of_octaves must be a positive integer.")
    if no_of_jobs == 1:
        power_trispectrum = get_single_observation_of_power_trispectrum(observable_string, output_directory,
                                                                        temperature_directory, beta, no_of_sites,
                                                                        no_of_equilibration_sweeps, no_of_octaves,
                                                                        base_time_period_shift, sampling_frequency)
    else:
        power_trispectra = pool.starmap(get_single_observation_of_power_trispectrum,
                                        [(observable_string, f"{output_directory}/job_{job_number + 1}",
                                          temperature_directory, beta, no_of_sites, no_of_equilibration_sweeps,
                                          no_of_octaves, base_time_period_shift, sampling_frequency)
                                         for job_number in range(no_of_jobs)])
        [print(index, power_trispectrum[0]) for index, power_trispectrum in enumerate(power_trispectra)]
        power_trispectrum = np.mean(np.array(power_trispectra, dtype=object), axis=0)
        print(power_trispectrum[0])
    # normalise estimator of power trispectrum with respect to its low-frequency value
    power_trispectrum[2] = [spectrum / spectrum[0] for spectrum in power_trispectrum[2]]
    return power_trispectrum


def get_power_trispectrum_direct(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                                 no_of_equilibration_sweeps, no_of_jobs, pool, no_of_octaves=2,
                                 base_time_period_shift=1, sampling_frequency=None):
    if no_of_octaves <= 0:
        raise Exception("no_of_octaves must be a positive integer.")
    if no_of_jobs == 1:
        sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
        power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(observable_string,
                                                                                    output_directory,
                                                                                    temperature_directory, beta,
                                                                                    no_of_sites,
                                                                                    no_of_equilibration_sweeps,
                                                                                    base_time_period_shift,
                                                                                    no_of_octaves, sampling_frequency)
    else:
        power_spectra_of_correlators = pool.starmap(get_power_spectra_of_trispectrum_correlators,
                                                    [(observable_string, f"{output_directory}/job_{job_number + 1}",
                                                      temperature_directory, beta, no_of_sites,
                                                      no_of_equilibration_sweeps, base_time_period_shift, no_of_octaves,
                                                      sampling_frequency)
                                                     for job_number in range(no_of_jobs)])
        power_spectra_of_correlators = np.mean(np.array(power_spectra_of_correlators, dtype=object), axis=0)
        sampling_frequency = get_sampling_frequency(f"{output_directory}/job_1", sampling_frequency,
                                                    temperature_directory)
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    spectra_in_frequency_shift_space = [np.absolute(item) for index, item in
                                        enumerate(np.fft.fft(transposed_power_spectra).transpose())
                                        if (index == 0 or index == 2 ** (math.floor(math.log(index, 2))))]
    # normalise power trispectrum with respect to its low-frequency value
    spectra_in_frequency_shift_space = [spectrum / spectrum[0] for spectrum in spectra_in_frequency_shift_space]
    return [np.atleast_1d(sampling_frequency / (2 ** (no_of_octaves + 1) * base_time_period_shift)),
            power_spectra_of_correlators[0, 0], spectra_in_frequency_shift_space]


def get_sampling_frequency(output_directory, sampling_frequency, temperature_directory):
    if sampling_frequency is None:
        acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        physical_time_scale = acceptance_rates[1] * acceptance_rates[0] ** 2 / 24.0
        sampling_frequency = 1.0 / physical_time_scale
    return sampling_frequency


def get_time_series(observable_string, output_directory, temperature_directory, beta, no_of_sites,
                    no_of_equilibration_sweeps):
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
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


def get_power_spectra_of_trispectrum_correlators(observable_string, output_directory, temperature_directory, beta,
                                                 no_of_sites, no_of_equilibration_sweeps, base_time_period_shift,
                                                 no_of_octaves, sampling_frequency):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    time_series = get_time_series(observable_string, output_directory, temperature_directory, beta, no_of_sites,
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


def get_single_observation_of_power_trispectrum(observable_string, output_directory, temperature_directory, beta,
                                                no_of_sites, no_of_equilibration_sweeps, no_of_octaves=2,
                                                base_time_period_shift=1, sampling_frequency=None):
    sampling_frequency = get_sampling_frequency(output_directory, sampling_frequency, temperature_directory)
    power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(observable_string, output_directory,
                                                                                temperature_directory, beta,
                                                                                no_of_sites, no_of_equilibration_sweeps,
                                                                                base_time_period_shift, no_of_octaves,
                                                                                sampling_frequency)
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    norm_of_spectra_in_frequency_shift_space = np.array(
        [np.absolute(item) for index, item in enumerate(np.fft.fft(transposed_power_spectra).transpose())
         if (index == 0 or index == 2 ** (math.floor(math.log(index, 2))))])
    return [np.atleast_1d(sampling_frequency / (2 ** (no_of_octaves + 1) * base_time_period_shift)),
            power_spectra_of_correlators[0, 0], norm_of_spectra_in_frequency_shift_space]
