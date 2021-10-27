from scipy import signal
import importlib
import math
import numpy as np

# import sample_getter module
sample_getter = importlib.import_module("sample_getter")


# main methods (see further down for a single-observation methods section and a basic methods section)

def get_power_spectrum(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                       no_of_equilibration_sweeps, no_of_jobs, pool, sampling_frequency=None):
    try:
        with open(f"{output_directory}/{observable_string}_power_spectrum_temp_eq_{temperature:.2f}.csv",
                  "r") as data_file:
            return np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        if no_of_jobs == 1:
            power_spectrum = get_single_observation_of_power_spectrum(algorithm_name, observable_string,
                                                                      output_directory, temperature, no_of_sites,
                                                                      no_of_equilibration_sweeps, sampling_frequency)
        else:
            power_spectra = pool.starmap(get_single_observation_of_power_spectrum,
                                         [(algorithm_name, observable_string,
                                           f"{output_directory}/job_{job_number + 1}", temperature, no_of_sites,
                                           no_of_equilibration_sweeps)
                                          for job_number in range(no_of_jobs)])
            power_spectrum = np.mean(np.array(power_spectra), axis=0)
        with open(f"{output_directory}/{observable_string}_power_spectrum_temp_eq_{temperature:.2f}.csv",
                  "w") as data_file:
            np.savetxt(data_file, power_spectrum, delimiter=",")
        return power_spectrum


def get_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                                     no_of_equilibration_sweeps, no_of_jobs, pool, time_period_shift=10,
                                     sampling_frequency=None):
    try:
        with open(f"{output_directory}/{observable_string}_power_spectrum_of_correlator_time_shift_eq_"
                  f"{time_period_shift}_delta_t_temp_eq_{temperature:.2f}.csv", "r") as data_file:
            return np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        if no_of_jobs == 1:
            correlator_power_spectrum = get_single_observation_of_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, time_period_shift, sampling_frequency)
        else:
            correlator_power_spectra = pool.starmap(
                get_single_observation_of_power_spectrum_of_correlator,
                [(algorithm_name, observable_string, f"{output_directory}/job_{job_number + 1}", temperature,
                  no_of_sites, no_of_equilibration_sweeps, time_period_shift)
                 for job_number in range(no_of_jobs)])
            correlator_power_spectrum = np.mean(np.array(correlator_power_spectra), axis=0)
        with open(f"{output_directory}/{observable_string}_power_spectrum_of_correlator_time_shift_eq_"
                  f"{time_period_shift}_delta_t_temp_eq_{temperature:.2f}.csv", "w") as data_file:
            np.savetxt(data_file, correlator_power_spectrum, delimiter=",")
        return correlator_power_spectrum


def get_power_trispectrum(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                          no_of_equilibration_sweeps, no_of_jobs, pool, no_of_auxiliary_frequency_octaves=2,
                          base_time_period_shift=1, sampling_frequency=None):
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        power_trispectrum = []
        stored_spectra = []
        with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_auxiliary_frequencies.csv", "r") as data_file:
            power_trispectrum.append(np.atleast_1d(np.loadtxt(data_file, dtype=float, delimiter=",")))
        for index in range(no_of_auxiliary_frequency_octaves + 1):
            with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                      f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                      f"{temperature:.2f}_f_prime_eq_{2 ** index}_x_delta_f_prime.csv", "r") as data_file:
                data = np.loadtxt(data_file, dtype=float, delimiter=",")
                if index == 0:
                    power_trispectrum.append(data[0])
                stored_spectra.append(data[1])
        power_trispectrum.append(np.array(stored_spectra))
    except IOError:
        if no_of_jobs == 1:
            power_trispectrum = get_single_observation_of_power_trispectrum(algorithm_name, observable_string,
                                                                            output_directory, temperature, no_of_sites,
                                                                            no_of_equilibration_sweeps,
                                                                            no_of_auxiliary_frequency_octaves,
                                                                            base_time_period_shift, sampling_frequency)
        else:
            power_trispectra = pool.starmap(get_single_observation_of_power_trispectrum,
                                            [(algorithm_name, observable_string,
                                              f"{output_directory}/job_{job_number + 1}", temperature, no_of_sites,
                                              no_of_equilibration_sweeps, no_of_auxiliary_frequency_octaves,
                                              base_time_period_shift, sampling_frequency)
                                             for job_number in range(no_of_jobs)])
            power_trispectrum = np.mean(np.array(power_trispectra, dtype=object), axis=0)
        with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_auxiliary_frequencies.csv", "w") as data_file:
            np.savetxt(data_file, power_trispectrum[0], delimiter=",")
        for index in range(no_of_auxiliary_frequency_octaves + 1):
            with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                      f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                      f"{temperature:.2f}_f_prime_eq_{2 ** index}_x_delta_f_prime.csv", "w") as data_file:
                np.savetxt(data_file, np.array([power_trispectrum[1], power_trispectrum[2][index]]), delimiter=",")
    return power_trispectrum


def get_power_trispectrum_zero_mode(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                                    no_of_equilibration_sweeps, no_of_jobs, pool, no_of_auxiliary_frequency_octaves=2,
                                    base_time_period_shift=1, sampling_frequency=None):
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_f_prime_eq_zero.csv", "r") as data_file:
            return np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        if no_of_jobs == 1:
            power_trispectrum_zero_mode = get_single_observation_of_power_trispectrum_zero_mode(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, no_of_auxiliary_frequency_octaves, base_time_period_shift,
                sampling_frequency)
        else:
            power_trispectra_zero_modes = pool.starmap(get_single_observation_of_power_trispectrum_zero_mode,
                                                       [(algorithm_name, observable_string,
                                                         f"{output_directory}/job_{job_number + 1}",
                                                         temperature, no_of_sites, no_of_equilibration_sweeps,
                                                         no_of_auxiliary_frequency_octaves, base_time_period_shift,
                                                         sampling_frequency) for job_number in range(no_of_jobs)])
            power_trispectrum_zero_mode = np.mean(np.array(power_trispectra_zero_modes, dtype=object), axis=0)
        with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_f_prime_eq_zero.csv", "w") as data_file:
            np.savetxt(data_file, power_trispectrum_zero_mode, delimiter=",")
        return power_trispectrum_zero_mode


def get_power_trispectrum_as_defined(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                                     no_of_equilibration_sweeps, no_of_jobs, pool, no_of_auxiliary_frequency_octaves=2,
                                     base_time_period_shift=1, sampling_frequency=None):
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        power_trispectrum = []
        stored_spectra = []
        with open(f"{output_directory}/{observable_string}_power_trispectrum_as_defined_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_auxiliary_frequencies.csv", "r") as data_file:
            power_trispectrum.append(np.atleast_1d(np.loadtxt(data_file, dtype=float, delimiter=",")))
        for index in range(no_of_auxiliary_frequency_octaves + 1):
            with open(f"{output_directory}/{observable_string}_power_trispectrum_as_defined_max_shift_eq_"
                      f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                      f"{temperature:.2f}_f_prime_eq_{2 ** index}_x_delta_f_prime.csv", "r") as data_file:
                data = np.loadtxt(data_file, dtype=float, delimiter=",")
                if index == 0:
                    power_trispectrum.append(data[0])
                stored_spectra.append(data[1])
        power_trispectrum.append(np.array(stored_spectra))
    except IOError:
        if no_of_jobs == 1:
            sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory,
                                                                      sampling_frequency, temperature)
            power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, base_time_period_shift, no_of_auxiliary_frequency_octaves,
                sampling_frequency)
        else:
            power_spectra_of_correlators = pool.starmap(get_power_spectra_of_trispectrum_correlators,
                                                        [(algorithm_name, observable_string,
                                                          f"{output_directory}/job_{job_number + 1}",
                                                          temperature, no_of_sites, no_of_equilibration_sweeps,
                                                          base_time_period_shift, no_of_auxiliary_frequency_octaves,
                                                          sampling_frequency) for job_number in range(no_of_jobs)])
            power_spectra_of_correlators = np.mean(np.array(power_spectra_of_correlators, dtype=object), axis=0)
            sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, f"{output_directory}/job_1",
                                                                      sampling_frequency, temperature)
        transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
        norm_of_spectra_in_auxiliary_frequency_space = [np.absolute(item) for index, item in
                                                        enumerate(np.fft.fft(transposed_power_spectra).transpose())
                                                        if 0 < index == 2 ** (math.floor(math.log(index, 2)))]
        power_trispectrum = [np.array([2 ** index for index in range(no_of_auxiliary_frequency_octaves + 1)])
                             * sampling_frequency / base_time_period_shift
                             / 2 ** (no_of_auxiliary_frequency_octaves + 1), power_spectra_of_correlators[0, 0],
                             norm_of_spectra_in_auxiliary_frequency_space]
        with open(f"{output_directory}/{observable_string}_power_trispectrum_as_defined_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_auxiliary_frequencies.csv", "w") as data_file:
            np.savetxt(data_file, power_trispectrum[0], delimiter=",")
        for index in range(no_of_auxiliary_frequency_octaves + 1):
            with open(f"{output_directory}/{observable_string}_power_trispectrum_as_defined_max_shift_eq_"
                      f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                      f"{temperature:.2f}_f_prime_eq_{2 ** index}_x_delta_f_prime.csv", "w") as data_file:
                np.savetxt(data_file, np.array([power_trispectrum[1], power_trispectrum[2][index]]), delimiter=",")
    return power_trispectrum


# single-observation methods

def get_single_observation_of_power_spectrum(algorithm_name, observable_string, output_directory, temperature,
                                             no_of_sites, no_of_equilibration_sweeps, sampling_frequency=None):
    sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory, sampling_frequency,
                                                              temperature)
    time_series = get_time_series(observable_string, output_directory, temperature, no_of_sites,
                                  no_of_equilibration_sweeps)
    return get_component_averaged_power_spectrum(time_series, sampling_frequency)


def get_single_observation_of_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory,
                                                           temperature, no_of_sites, no_of_equilibration_sweeps,
                                                           time_period_shift=10, sampling_frequency=None):
    sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory, sampling_frequency,
                                                              temperature)
    time_series = get_time_series(observable_string, output_directory, temperature, no_of_sites,
                                  no_of_equilibration_sweeps)
    if time_period_shift >= len(time_series[0]):
        raise Exception("time_period_shift must be an integer less than the sample size.")
    return get_component_averaged_power_spectrum(
        get_two_point_correlator(time_series - np.mean(time_series, axis=1), time_period_shift), sampling_frequency)


def get_single_observation_of_power_trispectrum(algorithm_name, observable_string, output_directory, temperature,
                                                no_of_sites, no_of_equilibration_sweeps,
                                                no_of_auxiliary_frequency_octaves=2, base_time_period_shift=1,
                                                sampling_frequency=None):
    sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory, sampling_frequency,
                                                              temperature)
    power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(algorithm_name, observable_string,
                                                                                output_directory, temperature,
                                                                                no_of_sites, no_of_equilibration_sweeps,
                                                                                base_time_period_shift,
                                                                                no_of_auxiliary_frequency_octaves,
                                                                                sampling_frequency)
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    norm_of_spectra_in_auxiliary_frequency_space = np.array(
        [np.absolute(item) for index, item in enumerate(np.fft.fft(transposed_power_spectra).transpose())
         if 0 < index == 2 ** (math.floor(math.log(index, 2)))])
    return [np.array([2 ** index for index in range(no_of_auxiliary_frequency_octaves + 1)]) * sampling_frequency
            / base_time_period_shift / 2 ** (no_of_auxiliary_frequency_octaves + 1), power_spectra_of_correlators[0, 0],
            norm_of_spectra_in_auxiliary_frequency_space]


def get_single_observation_of_power_trispectrum_zero_mode(algorithm_name, observable_string, output_directory,
                                                          temperature, no_of_sites, no_of_equilibration_sweeps,
                                                          no_of_auxiliary_frequency_octaves=2, base_time_period_shift=1,
                                                          sampling_frequency=None):
    power_spectra_of_trispectrum_correlators = get_power_spectra_of_trispectrum_correlators(
        algorithm_name, observable_string, output_directory, temperature, no_of_sites, no_of_equilibration_sweeps,
        base_time_period_shift, no_of_auxiliary_frequency_octaves, sampling_frequency)
    return [np.mean(power_spectra_of_trispectrum_correlators[:, 0], axis=0),
            np.sum(power_spectra_of_trispectrum_correlators[:, 1], axis=0)]


# basic methods

def get_time_series(observable_string, output_directory, temperature, no_of_sites, no_of_equilibration_sweeps):
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    sample = get_sample_method(output_directory, temperature, no_of_sites)[no_of_equilibration_sweeps:]
    sample = np.atleast_2d(sample)
    if len(sample) > 1:
        sample = sample.transpose()
    return sample


def get_component_averaged_power_spectrum(time_series, sampling_frequency):
    power_spectra = np.atleast_2d([signal.periodogram(component_2, fs=sampling_frequency) for component_2 in
                                   [component_1 - np.mean(component_1) for component_1 in time_series]])
    # now average over Cartesian components of original sample, where the `[:, 1:]' removes the f = 0 value as...
    # ...this value is invalid for a finite-time signal
    return np.mean(power_spectra, axis=0)[:, 1:]


def get_two_point_correlator(time_series, time_period_shift):
    """As the time series is not periodic, we chop off the first / last time_period_shift elements of each copy of the
        time series (rather than using np.roll()).  Previously, we returned
        np.conj(time_series[:, time_period_shift:]) * time_series[:, :len(time_series[0]) - time_period_shift],
        but have since removed the np.conj() operation as we only consider real-valued signals."""
    return time_series[:, time_period_shift:] * time_series[:, :len(time_series[0]) - time_period_shift]


def get_power_spectra_of_trispectrum_correlators(algorithm_name, observable_string, output_directory, temperature,
                                                 no_of_sites, no_of_equilibration_sweeps, base_time_period_shift,
                                                 no_of_auxiliary_frequency_octaves, sampling_frequency):
    sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory, sampling_frequency,
                                                              temperature)
    time_series = get_time_series(observable_string, output_directory, temperature, no_of_sites,
                                  no_of_equilibration_sweeps)
    if base_time_period_shift * 2 ** (no_of_auxiliary_frequency_octaves + 1) >= len(time_series[0]):
        raise Exception("base_time_period_shift * 2 ** (no_of_auxiliary_frequency_octaves + 1) must be less than the "
                        "smallest power of 2 that is less than the sample size.")
    # create 2 ** (no_of_auxiliary_frequency_octaves + 1) two-point correlators; here, [:, :len(time_series[0]) -
    # (2 ** (no_of_auxiliary_frequency_octaves + 1) - 1) * base_time_period_shift] ensures that (each component of) all
    # correlators are the same length, which ensures that their power spectra have common frequency values
    correlators = [
        get_two_point_correlator(np.atleast_2d([component - np.mean(component) for component in time_series]),
                                 i * base_time_period_shift)[
            :, :len(time_series[0]) - (2 ** (no_of_auxiliary_frequency_octaves + 1) - 1) * base_time_period_shift]
        for i in range(2 ** (no_of_auxiliary_frequency_octaves + 1))]
    return np.array(
        [get_component_averaged_power_spectrum(correlator, sampling_frequency) for correlator in correlators])
