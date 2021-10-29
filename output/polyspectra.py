"""This module contains methods that compute the polyspectra of some sample observable."""
from scipy import signal
import importlib
import math
import numpy as np
sample_getter = importlib.import_module("sample_getter")


"""Main methods (see further down for separate sections of single-observation methods and basic methods)"""


def get_power_spectrum(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                       no_of_equilibration_sweeps, no_of_jobs, pool, sampling_frequency=None):
    r"""
    Returns an estimate of the power spectrum S_X(f) := lim_{T -> inf} {E[| \Delta \tilde{X}_T(f) | ** 2] / T} of the
    time series X(t) of observable_string, where T is the total simulation time, \Delta \tilde{X}_T(f) is the Fourier
    transform of the truncated signal mean-zero time series
    \Delta X_T(t) := {X(t) - E[X] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) is considered a single observation of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero time series is
    \Delta \tilde{X}_T(f_k) := sum_{n = 0}^{N − 1} \Delta X(t_n) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_X(f_k) := lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the sample size (of the time series).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single observation of the power
    spectrum of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    no_of_jobs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The power spectrum.  A two-dimensional numpy array of shape (2, T / \Delta t / 2 - 1) /
        (2, (T / \Delta t - 1) / 2) for T / \Delta t even / odd.  Each element is a float.  The first / second
        sub-array is the frequencies / values of the power spectrum.
    """
    try:
        # first, attempt to open a previously computed estimate of the spectrum, then...
        with open(f"{output_directory}/{observable_string}_power_spectrum_temp_eq_{temperature:.2f}.csv",
                  "r") as data_file:
            return np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        # ...if the file does not exists, compute the estimate of the spectrum
        if no_of_jobs == 1:
            power_spectrum = get_single_observation_of_power_spectrum(algorithm_name, observable_string,
                                                                      output_directory, temperature, no_of_sites,
                                                                      no_of_equilibration_sweeps, sampling_frequency)
        else:
            # no_of_jobs > 1, so use the multiprocessing pool to compute the estimate of the spectrum corresponding to
            # each repeated simulation, then...
            power_spectra = pool.starmap(get_single_observation_of_power_spectrum,
                                         [(algorithm_name, observable_string,
                                           f"{output_directory}/job_{job_number + 1}", temperature, no_of_sites,
                                           no_of_equilibration_sweeps)
                                          for job_number in range(no_of_jobs)])
            # ...average over the results
            power_spectrum = np.mean(np.array(power_spectra), axis=0)
        # finally, save the estimated spectrum to file
        with open(f"{output_directory}/{observable_string}_power_spectrum_temp_eq_{temperature:.2f}.csv",
                  "w") as data_file:
            np.savetxt(data_file, power_spectrum, delimiter=",")
        return power_spectrum


def get_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                                     no_of_equilibration_sweeps, no_of_jobs, pool, time_period_shift=10,
                                     sampling_frequency=None):
    r"""
    Returns an estimate of the power spectrum S_Y(f, s) := lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T}
    of the correlator Y(t, s) := X(t) * X(t - s), where X(t) is the time series of observable_string,
    s = time_period_shift * \Delta t, \Delta t is the physical sampling interval, T is the total simulation time,
    \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) is considered a single observation of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * \Delta t / N}, where t_{n + 1} = t_n + \Delta t for all
    n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the discrete frequency spectrum and
    N = T / \Delta t is the sample size (of the correlator).  The factor of (\Delta t) ** 2 in the definition of the
    discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single observation of the (correlator)
    power spectrum of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    no_of_jobs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    time_period_shift : int, optional
        The number of multiples of the physical sampling interval by which the time series is shifted in order to form
        the correlator.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The power spectrum of the correlator.  A two-dimensional numpy array of shape (2, T / \Delta t / 2 - 1) /
        (2, (T / \Delta t - 1) / 2) for T / \Delta t even / odd.  Each element is a float.  The first / second
        sub-array is the frequencies / values of the power spectrum of the correlator.
    """
    try:
        # first, attempt to open a previously computed estimate of the spectrum, then...
        with open(f"{output_directory}/{observable_string}_power_spectrum_of_correlator_time_shift_eq_"
                  f"{time_period_shift}_delta_t_temp_eq_{temperature:.2f}.csv", "r") as data_file:
            return np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        # ...if the file does not exists, compute the estimate of the spectrum
        if no_of_jobs == 1:
            correlator_power_spectrum = get_single_observation_of_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, time_period_shift, sampling_frequency)
        else:
            # no_of_jobs > 1, so use the multiprocessing pool to compute the estimate of the spectrum corresponding to
            # each repeated simulation, then...
            correlator_power_spectra = pool.starmap(
                get_single_observation_of_power_spectrum_of_correlator,
                [(algorithm_name, observable_string, f"{output_directory}/job_{job_number + 1}", temperature,
                  no_of_sites, no_of_equilibration_sweeps, time_period_shift)
                 for job_number in range(no_of_jobs)])
            # ...average over the results
            correlator_power_spectrum = np.mean(np.array(correlator_power_spectra), axis=0)
        # finally, save the estimated spectrum to file
        with open(f"{output_directory}/{observable_string}_power_spectrum_of_correlator_time_shift_eq_"
                  f"{time_period_shift}_delta_t_temp_eq_{temperature:.2f}.csv", "w") as data_file:
            np.savetxt(data_file, correlator_power_spectrum, delimiter=",")
        return correlator_power_spectrum


def get_power_trispectrum(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                          no_of_equilibration_sweeps, no_of_jobs, pool, no_of_auxiliary_frequency_octaves=2,
                          base_time_period_shift=1, sampling_frequency=None):
    r"""
    Returns an estimate of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds of the
    time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total simulation
    time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency, f' is the
    auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a single observation of the
    dynamical 'distribution'.

    This shortcut estimator of the trispectrum in fact computes an estimate of
    E[|int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |], rather than of the
    direct definition above, which is encoded in get_power_trispectrum_as_defined().
    In conjunction with output/create_trispectrum_estimator_comparisons.py, the configuration file
    config_files/polyspectra/trispectrum_estimator_comparisons.txt shows that the current method is a low-noise
    equivalent of get_power_trispectrum_as_defined().

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the sample size (of the correlator).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single observation of the (correlator)
    power spectra of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    no_of_jobs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    !!!
        The power trispectrum.  !!!
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the trispectrum, then...
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
        # ...if the file does not exists, compute the estimate of the trispectrum
        if no_of_jobs == 1:
            power_trispectrum = get_single_observation_of_power_trispectrum(algorithm_name, observable_string,
                                                                            output_directory, temperature, no_of_sites,
                                                                            no_of_equilibration_sweeps,
                                                                            no_of_auxiliary_frequency_octaves,
                                                                            base_time_period_shift, sampling_frequency)
        else:
            # no_of_jobs > 1, so use the multiprocessing pool to compute the estimate of the trispectrum corresponding
            # to each repeated simulation, then...
            power_trispectra = pool.starmap(get_single_observation_of_power_trispectrum,
                                            [(algorithm_name, observable_string,
                                              f"{output_directory}/job_{job_number + 1}", temperature, no_of_sites,
                                              no_of_equilibration_sweeps, no_of_auxiliary_frequency_octaves,
                                              base_time_period_shift, sampling_frequency)
                                             for job_number in range(no_of_jobs)])
            # ...average over the results
            power_trispectrum = np.mean(np.array(power_trispectra, dtype=object), axis=0)
        # finally, save the estimated trispectrum to file
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
    r"""
    Returns an estimate of the zero (auxiliary-frequency) mode
    S_X^3(f, f' = 0) := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} ds of the power trispectrum of
    the time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total
    simulation time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency, f' is the
    auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a single observation of the
    dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the sample size (of the correlator).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single observation of the (correlator)
    power spectra of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    no_of_jobs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The zero mode of the power trispectrum.  A two-dimensional numpy array of shape (2, T / \Delta t / 2 - 1)
        / (2, (T / \Delta t - 1) / 2) for T / \Delta t even / odd.  Each element is a float.  The first / second
        sub-array is the frequencies / values of the zero mode of the power trispectrum.
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the trispectrum zero mode, then...
        with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_f_prime_eq_zero.csv", "r") as data_file:
            return np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        # ...if the file does not exists, compute the estimate of the trispectrum zero mode
        if no_of_jobs == 1:
            power_trispectrum_zero_mode = get_single_observation_of_power_trispectrum_zero_mode(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, no_of_auxiliary_frequency_octaves, base_time_period_shift,
                sampling_frequency)
        else:
            # no_of_jobs > 1, so use the multiprocessing pool to compute the estimate of the trispectrum zero mode
            # corresponding to each repeated simulation, then...
            power_trispectra_zero_modes = pool.starmap(get_single_observation_of_power_trispectrum_zero_mode,
                                                       [(algorithm_name, observable_string,
                                                         f"{output_directory}/job_{job_number + 1}",
                                                         temperature, no_of_sites, no_of_equilibration_sweeps,
                                                         no_of_auxiliary_frequency_octaves, base_time_period_shift,
                                                         sampling_frequency) for job_number in range(no_of_jobs)])
            # ...average over the results
            power_trispectrum_zero_mode = np.mean(np.array(power_trispectra_zero_modes, dtype=object), axis=0)
        # finally, save the estimated trispectrum zero mode to file
        with open(f"{output_directory}/{observable_string}_power_trispectrum_max_shift_eq_"
                  f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t_temp_eq_"
                  f"{temperature:.2f}_f_prime_eq_zero.csv", "w") as data_file:
            np.savetxt(data_file, power_trispectrum_zero_mode, delimiter=",")
        return power_trispectrum_zero_mode


def get_power_trispectrum_as_defined(algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                                     no_of_equilibration_sweeps, no_of_jobs, pool, no_of_auxiliary_frequency_octaves=2,
                                     base_time_period_shift=1, sampling_frequency=None):
    r"""
    Returns an estimate of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds of the
    time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total simulation
    time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency, f' is the
    auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a single observation of the
    dynamical 'distribution'.

    get_power_trispectrum() is a shortcut estimator of the trispectrum.  It computes an estimate of
    E[|int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |], rather than of the
    direct definition above, which is encoded in the current method.  In conjunction with
    output/create_trispectrum_estimator_comparisons.py, the configuration file
    config_files/polyspectra/trispectrum_estimator_comparisons.txt shows that get_power_trispectrum() is a low-noise
    equivalent of the current method.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the sample size (of the correlator).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single observation of the (correlator)
    power spectra of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    no_of_jobs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    !!!
        The power trispectrum.  !!!
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the trispectrum, then...
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
        # ...if the file does not exists, compute the estimate of the trispectrum
        if no_of_jobs == 1:
            sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory,
                                                                      sampling_frequency, temperature)
            power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, base_time_period_shift, no_of_auxiliary_frequency_octaves,
                sampling_frequency)
        else:
            # no_of_jobs > 1, so use the multiprocessing pool to compute the set of estimates of the spectra of the
            # correlators corresponding to each repeated simulation, then...
            power_spectra_of_correlators = pool.starmap(get_power_spectra_of_trispectrum_correlators,
                                                        [(algorithm_name, observable_string,
                                                          f"{output_directory}/job_{job_number + 1}",
                                                          temperature, no_of_sites, no_of_equilibration_sweeps,
                                                          base_time_period_shift, no_of_auxiliary_frequency_octaves,
                                                          sampling_frequency) for job_number in range(no_of_jobs)])
            # ...average over the set
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
        # finally, save the estimated trispectrum to file
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


"""Single-observation methods"""


def get_single_observation_of_power_spectrum(algorithm_name, observable_string, output_directory, temperature,
                                             no_of_sites, no_of_equilibration_sweeps, sampling_frequency=None):
    r"""
    Returns an estimate of a single observation lim_{T -> inf} [| \Delta \tilde{X}_T(f) | ** 2 / T] of the power
    spectrum S_X(f) := lim_{T -> inf} {E[| \Delta \tilde{X}_T(f) | ** 2] / T} of the time series X(t) of
    observable_string, where T is the total simulation time, \Delta \tilde{X}_T(f) is the Fourier transform of the
    truncated signal mean-zero time series \Delta X_T(t) := {X(t) - E[X] for all |t| <= T / 2, 0 otherwise}, t is time,
    f is frequency and E[.] is the expected value of the argument.  X(t) is considered a single observation of the
    dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero time series is
    \Delta \tilde{X}_T(f_k) := sum_{n = 0}^{N − 1} \Delta X(t_n) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_X(f_k) := lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the sample size (of the time series).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, the single observation of the power
    spectrum of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The single observation of the power spectrum.  A two-dimensional numpy array of shape
        (2, T / \Delta t / 2 - 1) / (2, (T / \Delta t - 1) / 2) for T / \Delta t even / odd.  Each element is a float.
        The first / second sub-array is the frequencies / values of the single observation of the power spectrum.
    """
    sampling_frequency = sample_getter.get_sampling_frequency(algorithm_name, output_directory, sampling_frequency,
                                                              temperature)
    time_series = get_time_series(observable_string, output_directory, temperature, no_of_sites,
                                  no_of_equilibration_sweeps)
    return get_component_averaged_power_spectrum(time_series, sampling_frequency)


def get_single_observation_of_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory,
                                                           temperature, no_of_sites, no_of_equilibration_sweeps,
                                                           time_period_shift=10, sampling_frequency=None):
    r"""
    Returns an estimate of a single observation lim_{T -> inf} [| \Delta \tilde{Y}_T(f) | ** 2 / T] of the power
    spectrum S_Y(f, s) := lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T}
    of the correlator Y(t, s) := X(t) * X(t - s), where X(t) is the time series of observable_string,
    s = time_period_shift * \Delta t, \Delta t is the physical sampling interval, T is the total simulation time,
    \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) is considered a single observation of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * \Delta t / N}, where t_{n + 1} = t_n + \Delta t for all
    n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the discrete frequency spectrum and
    N = T / \Delta t is the sample size (of the correlator).  The factor of (\Delta t) ** 2 in the definition of the
    discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, the single observation of the (correlator)
    power spectrum of each Cartesian component is computed and the average of the resultant quantities are returned.

    Parameters
    ----------
    algorithm_name : str
        The name of the sampling algorithm used to generate the time series / sample.
    observable_string : str
        The name of the observable whose power spectrum is to be estimated.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.
    time_period_shift : int, optional
        The number of multiples of the physical sampling interval by which the time series is shifted in order to form
        the correlator.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The single observation of the power spectrum of the correlator.  A two-dimensional numpy array of shape
        (2, T / \Delta t / 2 - 1) / (2, (T / \Delta t - 1) / 2) for T / \Delta t even / odd.  Each element is a float.
        The first / second sub-array is the frequencies / values of the single observation of the power spectrum of the
        correlator.
    """
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


"""Basic methods"""


def get_time_series(observable_string, output_directory, temperature, no_of_sites, no_of_equilibration_sweeps):
    r"""
    Returns the time series of observable_string.

    Parameters
    ----------
    observable_string : str
        The name of the observable whose time series is to be retrieved.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.

    Returns
    -------
    numpy.ndarray
        The time series / sample.  A two-dimensional numpy array of shape (n, T / \Delta t), where n >= 1 (an integer)
        is the number of Cartesian components of observable_string, T is the total simulation time and \Delta t is the
        sampling interval.
    """
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    sample = get_sample_method(output_directory, temperature, no_of_sites)[no_of_equilibration_sweeps:]
    sample = np.atleast_2d(sample)
    if len(sample) > 1:
        sample = sample.transpose()
    return sample


def get_component_averaged_power_spectrum(time_series, sampling_frequency):
    r"""
    Returns an estimate of the (Cartesian-)component average of a single observation
    lim_{T -> inf} [| \Delta \tilde{X}_T(f) | ** 2 / T] of the power spectrum
    S_X(f) := lim_{T -> inf} {E[| \Delta \tilde{X}_T(f) | ** 2] / T} of the time series X(t), where T is the total
    simulation time, \Delta \tilde{X}_T(f) is the Fourier transform of the truncated signal mean-zero time series
    \Delta X_T(t) := {X(t) - E[X] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) (a single component of time_series) is considered a single observation of the
    dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero time series is
    \Delta \tilde{X}_T(f_k) := sum_{n = 0}^{N − 1} \Delta X(t_n) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_X(f_k) := lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the sample size (of the time series).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    Parameters
    ----------
    time_series : numpy.ndarray
        The time series / sample whose component-averaged single-observation power spectrum is to be computed.  A
        two-dimensional numpy array of shape (n, T / \Delta t), where n >= 1 (an integer) is the number of components
        of time_series.
    sampling_frequency : float
        The sampling frequency.

    Returns
    -------
    numpy.ndarray
        The component average of the single observation of the power spectrum.  A two-dimensional numpy array of shape
        (2, T / \Delta t / 2 - 1) / (2, (T / \Delta t - 1) / 2) for T / \Delta t even / odd.  Each element is a float.
        The first / second sub-array is the frequencies / values of the component average of the single observation of
        the power spectrum.
    """
    power_spectra = np.atleast_2d([signal.periodogram(component_2, fs=sampling_frequency) for component_2 in
                                   [component_1 - np.mean(component_1) for component_1 in time_series]])
    # the `[:, 1:]' below removes the f = 0 value as this value is invalid for a finite-time stationary signal
    return np.mean(power_spectra, axis=0)[:, 1:]


def get_two_point_correlator(time_series, time_period_shift):
    r"""
    Returns the two-point correlator Y(t, s) := X(t) * X(t - s) of time_series X(t), where
    s = time_period_shift * \Delta t and \Delta t is the sampling interval.

    As the time series is not periodic, we chop off the first / last time_period_shift elements of each copy of the
    time series (rather than using np.roll()).  Previously, we returned
    np.conj(time_series[:, time_period_shift:]) * time_series[:, :len(time_series[0]) - time_period_shift],
    but have since removed the np.conj() operation as we only consider real-valued signals.

    Parameters
    ----------
    time_series : numpy.ndarray
        The time series / sample whose two-point correlator is to be computed.  A two-dimensional numpy array of shape
        (n, T / \Delta t), where n >= 1 (an integer) is the number of components of time_series.
    time_period_shift : int, optional
        The number of multiples of the physical sampling interval by which the time series is shifted in order to form
        the correlator.

    Returns
    -------
    numpy.ndarray
        The two-point correlator.  A two-dimensional numpy array of shape (n, T / \Delta t - 2 * time_period_shift),
        where n >= 1 (an integer) is the number of components of time_series and T is the total simulation time.
    """
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
