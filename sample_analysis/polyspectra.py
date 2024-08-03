"""This module contains methods that compute the polyspectra of some sample observable."""
from sample_getter import get_sampling_frequency
from scipy import signal
from typing import List
import math
import numpy as np
import sample_getter
import statistics


"""Main methods (see further down for separate sections of single-sample methods and base methods)"""


def get_power_spectrum(algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                       no_of_sites_string, external_global_moves_string, no_of_runs, pool,
                       no_of_equilibration_sweeps=None, thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of the power spectrum S_X(f) := lim_{T -> inf} {E[| \Delta \tilde{X}_T(f) | ** 2] / T} (with a
    standard error at each f) of the time series X(t) of observable_string, where T is the total simulation time,
    \Delta \tilde{X}_T(f) is the Fourier transform of the truncated signal mean-zero time series
    \Delta X_T(t) := {X(t) - E[X] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero time series is
    \Delta \tilde{X}_T(f_k) := sum_{n = 0}^{N − 1} \Delta X(t_n) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_X(f_k) := lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series.  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the power
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The power spectrum.  A two-dimensional numpy array of shape (3, N / 2 - 1) [(3, (N - 1) / 2)] for N even [odd].
        Each element is a float.  The first / second / third sub-array is the frequencies / values / standard errors of
        the power spectrum.  If N is even, the frequency spectrum is reduced to
        f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to
        f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the power spectrum is symmetric
        about f = 0 and its f = 0 value is invalid for a finite-time stationary signal.  If no_of_runs is 1, a numpy
        array of 1.0 floats is returned (as padding) for the standard errors.
    """
    try:
        # first, attempt to open a previously computed estimate of the spectrum...
        return np.load(f"{output_directory}/{observable_string}_power_spectrum_{algorithm_name.replace('-', '_')}_"
                       f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}.npy")
    except IOError:
        # ...then if the file does not exists, compute the estimate of the spectrum
        if no_of_runs == 1:
            power_spectrum = get_single_sample_of_power_spectrum(algorithm_name, observable_string,
                                                                      output_directory, temperature, temperature_index,
                                                                      no_of_sites, no_of_equilibration_sweeps,
                                                                      thinning_level, sampling_frequency)
            # append np.ones() to create the same shape array as for no_of_runs > 1 (where errors are the 3rd element)
            power_spectrum = np.concatenate([power_spectrum.flatten(),
                                             np.ones(len(power_spectrum[0]))]).reshape((3, -1))
        else:
            # no_of_runs > 1, so use the multiprocessing pool to compute the estimate of the spectrum corresponding to
            # each repeated simulation...
            power_spectra = pool.starmap(get_single_sample_of_power_spectrum,
                                         [(algorithm_name, observable_string,
                                           f"{output_directory}/run_{run_number}", temperature, temperature_index,
                                           no_of_sites, no_of_equilibration_sweeps, thinning_level)
                                          for run_number in range(no_of_runs)])
            # ...then compute the means and estimate the standard errors (wrt to the repeated simulations)
            power_spectrum = get_power_spectrum_means_and_std_errors(power_spectra)
        # finally, save the estimated spectrum to file
        np.save(f"{output_directory}/{observable_string}_power_spectrum_{algorithm_name.replace('-', '_')}_"
                f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}.npy", power_spectrum)
        return power_spectrum


def get_second_spectrum(algorithm_name, observable_string, output_directory, temperature, temperature_index,
                        no_of_sites, no_of_sites_string, external_global_moves_string, no_of_runs, pool,
                        no_of_equilibration_sweeps=None, thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of the second spectrum (minus its Gaussian contribution)
    S_Y(f, s = 0) := lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s = 0) | ** 2] / T} (with a standard error at each f) of
    the time series X(t) of observable_string, where Y(t, s = 0) := X(t) * X(t) is the square of the signal or the
    s = 0 correlator, T is the total simulation time, \Delta \tilde{Y}_T(f, s = 0) is the Fourier transform of the
    truncated mean-zero square of the signal
    \Delta Y_T(t, s = 0) := {Y(t, s = 0) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.]
    is the expected value of the argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero square of the signal is
    \tilde{Y}_T(f_k, s = 0) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s = 0) exp(- 2 * pi * i * f_k * t_n), so that the
    estimate of the second spectrum is
    S_Y(f_k, s = 0) := lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s = 0) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s = 0) ∣ ** 2] * \Delta t / N}, where t_{n + 1} = t_n + \Delta t for
    all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the discrete frequency spectrum and
    N = T / \Delta t is the number of samples of the time series (of the second spectrum).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the second
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The second spectrum.  A two-dimensional numpy array of shape (3, N / 2 - 1) [(3, (N - 1) / 2)] for N even
        [odd].  Each element is a float.  The first / second / third sub-array is the frequencies / values /
        standard errors of the second spectrum.  If N is even, the frequency spectrum is reduced to
        f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to
        f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the power spectrum is symmetric
        about f = 0 and its f = 0 value is invalid for a finite-time stationary signal.  If no_of_runs is 1, a numpy
        array of 1.0 floats is returned (as padding) for the standard errors.
    """
    return get_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory, temperature,
                                            temperature_index, no_of_sites, no_of_sites_string,
                                            external_global_moves_string, no_of_runs, pool, 0,
                                            no_of_equilibration_sweeps, thinning_level, sampling_frequency)


def get_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory, temperature,
                                     temperature_index, no_of_sites, no_of_sites_string, external_global_moves_string,
                                     no_of_runs, pool, time_period_shift=10, no_of_equilibration_sweeps=None,
                                     thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of the power spectrum S_Y(f, s) := lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T}
    (with a standard error at each f) of the correlator Y(t, s) := X(t) * X(t - s), where X(t) is the time series of
    observable_string, s = time_period_shift * \Delta t, \Delta t is the physical sampling interval, T is the total
    simulation time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * \Delta t / N}, where t_{n + 1} = t_n + \Delta t for all
    n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the discrete frequency spectrum and
    N = T / \Delta t is the number of samples of the time series (of the correlator).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    time_period_shift : int, optional
        The number of multiples of the physical sampling interval by which the time series is shifted in order to form
        the correlator.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The power spectrum of the correlator.  A two-dimensional numpy array of shape (3, N / 2 - 1) [(3, (N - 1) / 2)]
        for N even [odd].  Each element is a float.  The first / second / third sub-array is the frequencies / values /
        standard errors of the power spectrum of the correlator.  If N is even, the frequency spectrum is reduced to
        f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to
        f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the power spectrum is symmetric
        about f = 0 and its f = 0 value is invalid for a finite-time stationary signal.  If no_of_runs is 1, a numpy
        array of 1.0 floats is returned (as padding) for the standard errors.
    """
    try:
        # first, attempt to open a previously computed estimate of the spectrum...
        return np.load(f"{output_directory}/{observable_string}_power_spectrum_of_correlator_"
                       f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq"
                       f"_{temperature:.4f}_time_shift_eq_{time_period_shift}_delta_t.npy")
    except IOError:
        # ...then if the file does not exists, compute the estimate of the spectrum
        if no_of_runs == 1:
            correlator_power_spectrum = get_single_sample_of_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                time_period_shift, no_of_equilibration_sweeps, thinning_level, sampling_frequency)
            # append np.ones() to create the same shape array as for no_of_runs > 1 (where errors are the 3rd element)
            correlator_power_spectrum = np.concatenate([correlator_power_spectrum.flatten(),
                                                        np.ones(len(correlator_power_spectrum[0]))]).reshape((3, -1))
        else:
            # no_of_runs > 1, so use the multiprocessing pool to compute the estimate of the spectrum corresponding to
            # each repeated simulation...
            correlator_power_spectra = pool.starmap(
                get_single_sample_of_power_spectrum_of_correlator,
                [(algorithm_name, observable_string, f"{output_directory}/run_{run_number}", temperature,
                  temperature_index, no_of_sites, time_period_shift, no_of_equilibration_sweeps, thinning_level)
                 for run_number in range(no_of_runs)])
            # ...then compute the means and estimate the standard errors (wrt to the repeated simulations)
            correlator_power_spectrum = get_power_spectrum_means_and_std_errors(correlator_power_spectra)
        # finally, save the estimated spectrum to file
        np.save(f"{output_directory}/{observable_string}_power_spectrum_of_correlator_"
                f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
                f"{temperature:.4f}_time_shift_eq_{time_period_shift}_delta_t.npy", correlator_power_spectrum)
        return correlator_power_spectrum


def get_power_trispectrum(algorithm_name, observable_string, output_directory, temperature, temperature_index,
                          no_of_sites, no_of_sites_string, external_global_moves_string, no_of_runs, pool,
                          no_of_auxiliary_frequency_octaves=6, base_time_period_shift=1,
                          no_of_equilibration_sweeps=None, thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of the shortcut estimator
    E[| int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |] of the complex norm
    of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds (with
    standard errors at each f) of the time series X(t) of observable_string, where the correlator
    Y(t, s) := X(t) * X(t - s), T is the total simulation time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of
    the truncated mean-zero correlator \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is
    time, s is the auxiliary time, f is frequency, f' is the auxiliary frequency and E[.] is the expected value of the
    argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    In conjunction with output/create_trispectrum_estimator_comparisons.py, the configuration file
    config_files/polyspectra/trispectrum_estimator_comparisons.txt shows that this shortcut estimator (the current
    method) is a low-noise equivalent of the direct definition (of the complex norm) encoded in
    get_power_trispectrum_as_defined().

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    List[numpy.ndarray]
        The power trispectrum.  A list of length 4.  The first component is the auxiliary frequencies and is a
        one-dimensional numpy (of floats) of length no_of_auxiliary_frequency_octaves + 1.  The second component is the
        frequencies and is a one-dimensional numpy array (of floats) of length N / 2 - 1 [(N - 1) / 2] for N even
        [odd].  If N is even, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the correlator power spectra are
        symmetric about f = 0 and their f = 0 values are invalid for a finite-time stationary signal.  The third
        component is a two-dimensional numpy array (of floats) of shape
        (no_of_auxiliary_frequency_octaves + 1, N / 2 - 1) [(no_of_auxiliary_frequency_octaves + 1, (N - 1) / 2)] for N
        even [odd].  The nth sub-array of the third / fourth component is the trispectrum / trispectrum standard errors
        at the auxiliary-frequency value given by the nth element of the first component.  If no_of_runs is 1, a numpy
        array of 1.0 floats is returned (as padding) for the standard errors.
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the trispectrum...
        auxiliary_frequencies = np.loadtxt(
            f"{output_directory}/{observable_string}_trispectrum_auxiliary_frequencies_"
            f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
            f"{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_"
            f"delta_t.csv", dtype=float, delimiter=",")
        npz_file = np.load(f"{output_directory}/{observable_string}_trispectrum_{algorithm_name.replace('-', '_')}_"
                           f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq"
                           f"_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npz")
        return [auxiliary_frequencies, npz_file["arr_0"][0],
                np.array([npz_file[f"arr_{index}"][1] for index in range(no_of_auxiliary_frequency_octaves + 1)]),
                np.array([npz_file[f"arr_{index}"][2] for index in range(no_of_auxiliary_frequency_octaves + 1)])]
    except IOError:
        # ...then if the file does not exists, compute the estimate of the trispectrum
        if no_of_runs == 1:
            power_trispectrum = get_single_sample_of_power_trispectrum(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                no_of_auxiliary_frequency_octaves, base_time_period_shift, no_of_equilibration_sweeps, thinning_level,
                sampling_frequency)
            # append np.ones() to create the same shape array as for no_of_runs > 1 (where errors are the 4th element)
            power_trispectrum.append(np.ones(np.shape(power_trispectrum[2])))
        else:
            # no_of_runs > 1, so use the multiprocessing pool to compute the estimate of the trispectrum corresponding
            # to each repeated simulation...
            power_trispectra = pool.starmap(get_single_sample_of_power_trispectrum,
                                            [(algorithm_name, observable_string,
                                              f"{output_directory}/run_{run_number}", temperature,
                                              temperature_index, no_of_sites, no_of_auxiliary_frequency_octaves,
                                              base_time_period_shift, no_of_equilibration_sweeps, thinning_level,
                                              sampling_frequency) for run_number in range(no_of_runs)])
            # ...then extract the frequencies and spectra...
            auxiliary_frequencies = np.mean([simulation[0] for simulation in power_trispectra], axis=0)
            frequencies = np.mean(np.array([simulation[1] for simulation in power_trispectra]), axis=0)
            trispectra_sans_any_frequencies = np.array([simulation[2] for simulation in power_trispectra])
            # ...then average over results and estimate standard errors
            power_trispectrum = [auxiliary_frequencies, frequencies, np.mean(trispectra_sans_any_frequencies, axis=0),
                                 np.std(trispectra_sans_any_frequencies, axis=0) / len(power_trispectra) ** 0.5]
        # finally, save the estimated trispectrum to file
        np.savetxt(f"{output_directory}/{observable_string}_trispectrum_auxiliary_frequencies_"
                   f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
                   f"{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}"
                   f"_delta_t.csv", power_trispectrum[0], delimiter=",")
        np.savez(f"{output_directory}/{observable_string}_trispectrum_{algorithm_name.replace('-', '_')}_"
                 f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq_"
                 f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npz",
                 *[np.array([power_trispectrum[1], power_trispectrum[2][index], power_trispectrum[3][index]])
                   for index in range(no_of_auxiliary_frequency_octaves + 1)])
        return power_trispectrum


def get_power_trispectrum_zero_mode(algorithm_name, observable_string, output_directory, temperature, temperature_index,
                                    no_of_sites, no_of_sites_string, external_global_moves_string, no_of_runs, pool,
                                    no_of_auxiliary_frequency_octaves=6, base_time_period_shift=1,
                                    no_of_equilibration_sweeps=None, thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of the zero (auxiliary-frequency) mode
    S_X^3(f, f' = 0) := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} ds of the power trispectrum (with
    standard errors at each f) of the time series X(t) of observable_string, where the correlator
    Y(t, s) := X(t) * X(t - s), T is the total simulation time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of
    the truncated mean-zero correlator \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is
    time, s is the auxiliary time, f is frequency, f' is the auxiliary frequency and E[.] is the expected value of the
    argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The zero mode of the power trispectrum.  A two-dimensional numpy array of shape (3, N / 2 - 1)
        [(2, (N - 1) / 2)] for N even [odd].  Each element is a float.  The first / second / third sub-array is the
        frequencies / values / standard errors of the zero mode of the power trispectrum.  If N is even, the frequency
        spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency
        spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the
        correlator power spectra are symmetric about f = 0 and their f = 0 values are invalid for a finite-time
        stationary signal.  If no_of_runs is 1, a numpy array of 1.0 floats is returned (as padding) for the standard
        errors.
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the trispectrum zero mode...
        return np.load(f"{output_directory}/{observable_string}_trispectrum_zero_mode_"
                       f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq"
                       f"_{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_"
                       f"{base_time_period_shift}_delta_t.npy")
    except IOError:
        # ...then if the file does not exists, compute the estimate of the trispectrum zero mode
        if no_of_runs == 1:
            power_trispectrum_zero_mode = get_single_sample_of_power_trispectrum_zero_mode(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                no_of_auxiliary_frequency_octaves, base_time_period_shift, no_of_equilibration_sweeps, thinning_level,
                sampling_frequency)
            # append np.ones() to create the same shape array as for no_of_runs > 1 (where errors are the 3rd element)
            power_trispectrum_zero_mode = np.concatenate(
                [power_trispectrum_zero_mode.flatten(), np.ones(len(power_trispectrum_zero_mode[0]))]).reshape((3, -1))
        else:
            # no_of_runs > 1, so use the multiprocessing pool to compute the estimate of the trispectrum zero mode
            # corresponding to each repeated simulation...
            power_trispectra_zero_modes = pool.starmap(get_single_sample_of_power_trispectrum_zero_mode,
                                                       [(algorithm_name, observable_string,
                                                         f"{output_directory}/run_{run_number}",
                                                         temperature, temperature_index, no_of_sites,
                                                         no_of_auxiliary_frequency_octaves, base_time_period_shift,
                                                         no_of_equilibration_sweeps, thinning_level, sampling_frequency)
                                                        for run_number in range(no_of_runs)])
            # ...then extract the frequencies and spectra...
            frequencies = np.mean(np.array([simulation[0] for simulation in power_trispectra_zero_modes]), axis=0)
            trispectra_sans_frequencies = np.array([simulation[1] for simulation in power_trispectra_zero_modes])
            # ...then average over results and estimate standard errors
            power_trispectrum_zero_mode = np.array(
                [frequencies, np.mean(trispectra_sans_frequencies, axis=0),
                 np.std(trispectra_sans_frequencies, axis=0) / len(power_trispectra_zero_modes) ** 0.5])
        # finally, save the estimated trispectrum zero mode to file
        np.save(f"{output_directory}/{observable_string}_trispectrum_zero_mode_{algorithm_name.replace('-', '_')}_"
                f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq_"
                f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npy",
                power_trispectrum_zero_mode)
        return power_trispectrum_zero_mode


def get_power_trispectrum_nonzero_mode(algorithm_name, observable_string, output_directory, temperature,
                                       temperature_index, no_of_sites, no_of_sites_string, external_global_moves_string,
                                       no_of_runs, pool, no_of_auxiliary_frequency_octaves=6, base_time_period_shift=1,
                                       target_auxiliary_frequency=None, no_of_equilibration_sweeps=None,
                                       thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of a single auxiliary-frequency mode of the shortcut estimator
    E[| int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |] of the complex norm
    of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds (with
    standard errors at each f) of the time series X(t) of observable_string, where the correlator
    Y(t, s) := X(t) * X(t - s), T is the total simulation time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of
    the truncated mean-zero correlator \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is
    time, s is the auxiliary time, f is frequency, f' is the auxiliary frequency and E[.] is the expected value of the
    argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    In conjunction with output/create_trispectrum_estimator_comparisons.py, the configuration file
    config_files/polyspectra/trispectrum_estimator_comparisons.txt shows that this shortcut estimator (the current
    method) is a low-noise equivalent of the direct definition (of the complex norm) encoded in
    get_power_trispectrum_as_defined().

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    target_auxiliary_frequency : float or None, optional
        The target auxiliary frequency.  Within each repeated simulation, the auxiliary-frequency mode (of the
        trispectrum) whose auxiliary-frequency value is closest to target_auxiliary_frequency is chosen.  If None is
        given, the largest auxiliary frequency is chosen.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    List[numpy.ndarray]
        Single mode of the power trispectrum.  A list of length 4.  The first component is the auxiliary frequency and
        is a float.  The second, third and fourth components are one-dimensional numpy arrays (of floats) of length
        N / 2 - 1 [(N - 1) / 2] for N even [odd].  The second component is the frequencies.  If N is even, the
        frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the
        frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because
        the correlator power spectra are symmetric about f = 0 and their f = 0 values are invalid for a finite-time
        stationary signal.  The third / fourth component is the trispectrum / trispectrum standard errors at the
        auxiliary-frequency value given by the first component.  If no_of_runs is 1, a numpy array of 1.0 floats is
        returned (as padding) for the standard errors.
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the chosen mode of the trispectrum...
        auxiliary_frequency = np.loadtxt(
            f"{output_directory}/{observable_string}_trispectrum_nonzero_mode_auxiliary_frequency_"
            f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
            f"{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_"
            f"delta_t.csv", dtype=float, delimiter=",")
        frequencies, spectrum, errors = np.load(
            f"{output_directory}/{observable_string}_trispectrum_nonzero_mode_{algorithm_name.replace('-', '_')}_"
            f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq_"
            f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npy")
        return [auxiliary_frequency, frequencies, spectrum, errors]
    except IOError:
        # ...then if the files do not exist, compute the estimate of the nonzero mode of the trispectrum
        if no_of_runs == 1:
            power_trispectrum_nonzero_mode = get_single_sample_of_power_trispectrum_nonzero_mode(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                no_of_auxiliary_frequency_octaves, base_time_period_shift, target_auxiliary_frequency,
                no_of_equilibration_sweeps, thinning_level, sampling_frequency)
            # append np.ones() to create the same shape array as for no_of_runs > 1 (where errors are the 4th element)
            power_trispectrum_nonzero_mode.append(np.ones(np.shape(power_trispectrum_nonzero_mode[2])))
        else:
            # no_of_runs > 1, so use the multiprocessing pool to compute the estimate of the trispectrum corresponding
            # to each repeated simulation...
            power_trispectra_nonzero_mode = pool.starmap(
                get_single_sample_of_power_trispectrum_nonzero_mode,
                [(algorithm_name, observable_string, f"{output_directory}/run_{run_number}", temperature,
                  temperature_index, no_of_sites, no_of_auxiliary_frequency_octaves, base_time_period_shift,
                  target_auxiliary_frequency, no_of_equilibration_sweeps, thinning_level, sampling_frequency)
                 for run_number in range(no_of_runs)])
            # ...then extract the frequencies and spectra...
            auxiliary_frequency = statistics.mean([simulation[0] for simulation in power_trispectra_nonzero_mode])
            frequencies = np.mean(np.array([simulation[1] for simulation in power_trispectra_nonzero_mode]), axis=0)
            trispectra_sans_any_frequencies = np.array([simulation[2] for simulation in power_trispectra_nonzero_mode])
            # ...then average over results and estimate standard errors
            power_trispectrum_nonzero_mode = [
                auxiliary_frequency, frequencies, np.mean(trispectra_sans_any_frequencies, axis=0),
                np.std(trispectra_sans_any_frequencies, axis=0) / len(power_trispectra_nonzero_mode) ** 0.5]
        # finally, save the estimated nonzero mode of trispectrum to file
        np.savetxt(f"{output_directory}/{observable_string}_trispectrum_nonzero_mode_auxiliary_frequency_"
                   f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
                   f"{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}"
                   f"_delta_t.csv", np.atleast_1d(power_trispectrum_nonzero_mode[0]), delimiter=",")
        np.save(f"{output_directory}/{observable_string}_trispectrum_nonzero_mode_{algorithm_name.replace('-', '_')}_"
                f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq_"
                f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npy",
                np.array([power_trispectrum_nonzero_mode[1], power_trispectrum_nonzero_mode[2],
                          power_trispectrum_nonzero_mode[3]]))
        return power_trispectrum_nonzero_mode


def get_power_trispectrum_as_defined(algorithm_name, observable_string, output_directory, temperature,
                                     temperature_index, no_of_sites, no_of_sites_string, external_global_moves_string,
                                     no_of_runs, pool, no_of_auxiliary_frequency_octaves=6, base_time_period_shift=1,
                                     no_of_equilibration_sweeps=None, thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of the complex norm (as defined) of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds of the
    time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total simulation
    time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, s is the auxiliary time, f is
    frequency, f' is the auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a
    single sample of the dynamical 'distribution'.

    get_power_trispectrum() encodes the shortcut estimator
    E[| int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |] of the complex norm
    of the trispectrum, rather than the direct definition, which is encoded in the current method.  In conjunction with
    output/create_trispectrum_estimator_comparisons.py, the configuration file
    config_files/polyspectra/trispectrum_estimator_comparisons.txt shows that get_power_trispectrum() is a low-noise
    equivalent of the current method.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_sites_string : str
        The string describing the number of lattice sites.
    external_global_moves_string : str
        A string that describes whether or not global topological-sector or twist moves were employed in Markov process.
    no_of_runs : int
        The number of repeated simulations.
    pool : multiprocessing.Pool() or None
        The multiprocessing pool (for multiple repeated simulations) or None (for a single simulation).
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    List[numpy.ndarray]
        The power trispectrum (as defined).  A list of length 3.  The first component is the auxiliary frequencies and
        is a one-dimensional numpy (of floats) of length no_of_auxiliary_frequency_octaves + 1.  The second component
        is the frequencies and is a one-dimensional numpy array (of floats) of length N / 2 - 1 [(N - 1) / 2] for N
        even [odd].  If N is even, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the correlator power spectra are
        symmetric about f = 0 and their f = 0 values are invalid for a finite-time stationary signal.  The third
        component is a two-dimensional numpy array (of floats) of shape
        (no_of_auxiliary_frequency_octaves + 1, N / 2 - 1) [(no_of_auxiliary_frequency_octaves + 1, (N - 1) / 2)] for N
        even [odd].  The nth sub-array of the third component is the trispectrum at the auxiliary-frequency value given
        by the nth element of the first component.
    """
    if no_of_auxiliary_frequency_octaves <= 0:
        raise Exception("no_of_auxiliary_frequency_octaves must be a positive integer.")
    try:
        # first, attempt to open a previously computed estimate of the trispectrum...
        auxiliary_frequencies = np.atleast_1d(np.loadtxt(
            f"{output_directory}/{observable_string}_trispectrum_as_defined_auxiliary_frequencies_"
            f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
            f"{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_"
            f"delta_t.csv", dtype=float, delimiter=","))
        npz_file = np.load(f"{output_directory}/{observable_string}_trispectrum_as_defined_"
                           f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_"
                           f"{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq_"
                           f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npz",
                           allow_pickle=True)
        return [auxiliary_frequencies, npz_file["arr_0"][0],
                np.array([npz_file[f"arr_{index}"][1] for index in range(no_of_auxiliary_frequency_octaves + 1)])]
    except IOError:
        # ...then if the file does not exists, compute the estimate of the trispectrum
        if no_of_runs == 1:
            sampling_frequency = get_sampling_frequency(algorithm_name, output_directory, sampling_frequency,
                                                        temperature_index)
            power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                base_time_period_shift, no_of_auxiliary_frequency_octaves, no_of_equilibration_sweeps, thinning_level,
                sampling_frequency)
        else:
            sampling_frequency = get_sampling_frequency(algorithm_name, f"{output_directory}/run_1", sampling_frequency,
                                                        temperature_index)
            # no_of_runs > 1, so use the multiprocessing pool to compute the set of estimates of the spectra of the
            # correlators corresponding to each repeated simulation...
            power_spectra_of_correlators = pool.starmap(get_power_spectra_of_trispectrum_correlators, [
                (algorithm_name, observable_string, f"{output_directory}/run_{run_number}", temperature,
                 temperature_index, no_of_sites, base_time_period_shift, no_of_auxiliary_frequency_octaves,
                 no_of_equilibration_sweeps, thinning_level, sampling_frequency) for run_number in range(no_of_runs)])
            # ...then average over the set
            power_spectra_of_correlators = np.mean(np.array(power_spectra_of_correlators, dtype=object), axis=0)
        transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
        norm_of_spectra_in_auxiliary_frequency_space = np.array(
            [np.absolute(item) for index, item in enumerate(np.fft.fft(transposed_power_spectra).transpose())
             if 0 < index == 2 ** (math.floor(math.log(index, 2)))])
        power_trispectrum = [np.array([2 ** index for index in range(no_of_auxiliary_frequency_octaves + 1)])
                             * sampling_frequency / base_time_period_shift
                             / 2 ** (no_of_auxiliary_frequency_octaves + 1), power_spectra_of_correlators[0, 0],
                             norm_of_spectra_in_auxiliary_frequency_space]
        # finally, save the estimated trispectrum to file
        np.savetxt(f"{output_directory}/{observable_string}_trispectrum_as_defined_auxiliary_frequencies_"
                   f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq_"
                   f"{temperature:.4f}_max_shift_eq_{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}"
                   f"_delta_t.csv", power_trispectrum[0], delimiter=",")
        np.savez(f"{output_directory}/{observable_string}_trispectrum_as_defined_{algorithm_name.replace('-', '_')}_"
                 f"{external_global_moves_string}_{no_of_sites_string}_temp_eq_{temperature:.4f}_max_shift_eq_"
                 f"{2 ** no_of_auxiliary_frequency_octaves}_x_{base_time_period_shift}_delta_t.npz",
                 *[np.array([power_trispectrum[1], power_trispectrum[2][index]])
                   for index in range(no_of_auxiliary_frequency_octaves + 1)])
        return power_trispectrum


"""Single-sample methods"""


def get_single_sample_of_power_spectrum(algorithm_name, observable_string, output_directory, temperature,
                                             temperature_index, no_of_sites, no_of_equilibration_sweeps=None,
                                             thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of a single sample lim_{T -> inf} [| \Delta \tilde{X}_T(f) | ** 2 / T] of the power
    spectrum S_X(f) := lim_{T -> inf} {E[| \Delta \tilde{X}_T(f) | ** 2] / T} of the time series X(t) of
    observable_string, where T is the total simulation time, \Delta \tilde{X}_T(f) is the Fourier transform of the
    truncated signal mean-zero time series \Delta X_T(t) := {X(t) - E[X] for all |t| <= T / 2, 0 otherwise}, t is time,
    f is frequency and E[.] is the expected value of the argument.  X(t) is considered a single sample of the
    dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero time series is
    \Delta \tilde{X}_T(f_k) := sum_{n = 0}^{N − 1} \Delta X(t_n) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_X(f_k) := lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series.  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, the single sample of the power
    spectrum of each Cartesian component is computed and the average of the resultant quantities is returned.

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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The single sample of the power spectrum.  A two-dimensional numpy array of shape (2, N / 2 - 1)
        [(2, (N - 1) / 2)] for N even [odd].  Each element is a float.  The first / second sub-array is the frequencies
        / values of the power spectrum.  If N is even, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the power spectrum is symmetric about
        f = 0 and its f = 0 value is invalid for a finite-time stationary signal.
    """
    sampling_frequency = get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_index)
    time_series = get_time_series(observable_string, output_directory, temperature, temperature_index, no_of_sites,
                                  no_of_equilibration_sweeps, thinning_level)
    return get_component_averaged_power_spectrum(time_series, sampling_frequency)


def get_single_sample_of_power_spectrum_of_correlator(algorithm_name, observable_string, output_directory,
                                                           temperature, temperature_index, no_of_sites,
                                                           time_period_shift=10, no_of_equilibration_sweeps=None,
                                                           thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of a single sample lim_{T -> inf} [| \Delta \tilde{Y}_T(f) | ** 2 / T] of the power
    spectrum S_Y(f, s) := lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T}
    of the correlator Y(t, s) := X(t) * X(t - s), where X(t) is the time series of observable_string,
    s = time_period_shift * \Delta t, \Delta t is the physical sampling interval, T is the total simulation time,
    \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{Y}_T(f_k, s) ∣ ** 2] * \Delta t / N}, where t_{n + 1} = t_n + \Delta t for all
    n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the discrete frequency spectrum and
    N = T / \Delta t is the number of samples of the time series (of the correlator).  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, the single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    time_period_shift : int, optional
        The number of multiples of the physical sampling interval by which the time series is shifted in order to form
        the correlator.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The single sample of the power spectrum of the correlator.  A two-dimensional numpy array of shape
        (2, N / 2 - 1) [(2, (N - 1) / 2)] for N even [odd].  Each element is a float.  The first / second sub-array is
        the frequencies / values of the single sample of the power spectrum of the correlator.  If N is even, the
        frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the
        frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because
        the power spectrum is symmetric about f = 0 and its f = 0 value is invalid for a finite-time stationary signal.
    """
    sampling_frequency = get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_index)
    time_series = get_time_series(observable_string, output_directory, temperature, temperature_index, no_of_sites,
                                  no_of_equilibration_sweeps, thinning_level)
    if time_period_shift >= len(time_series[0]):
        raise Exception("time_period_shift must be an integer less than the sample size.")
    return get_component_averaged_power_spectrum(get_two_point_correlator(
        np.atleast_2d([component - np.mean(component) for component in time_series]), time_period_shift),
        sampling_frequency)


def get_single_sample_of_power_trispectrum(algorithm_name, observable_string, output_directory, temperature,
                                                temperature_index, no_of_sites, no_of_auxiliary_frequency_octaves=6,
                                                base_time_period_shift=1, no_of_equilibration_sweeps=None,
                                                thinning_level=None, sampling_frequency=None):
    r"""
    Returns an estimate of a single sample
    | int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) | ds of the shortcut
    estimator E[| int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |] of the
    complex norm of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds of the
    time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total simulation
    time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, s is the auxiliary time, f is
    frequency, f' is the auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a
    single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    List[numpy.ndarray]
        The single sample of the power trispectrum.  A list of length 3.  The first component is the auxiliary
        frequencies and is a one-dimensional numpy array (of floats) of length no_of_auxiliary_frequency_octaves + 1.
        The second component is the frequencies and is a one-dimensional numpy array (of floats) of length N / 2 - 1
        [(N - 1) / 2] for N even [odd].  If N is even, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is reduced to f_k \in
        {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the correlator power spectra are
        symmetric about f = 0 and their f = 0 values are invalid for a finite-time stationary signal.  The third
        component is a two-dimensional numpy array (of floats) of shape
        (no_of_auxiliary_frequency_octaves + 1, N / 2 - 1) [(no_of_auxiliary_frequency_octaves + 1, N / 2)] for N even
        [odd].  The nth sub-array of the third component is the (single sample of the shortcut estimator of the)
        trispectrum at the auxiliary-frequency value given by the nth element of the first component.
    """
    sampling_frequency = get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_index)
    power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(
        algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
        base_time_period_shift, no_of_auxiliary_frequency_octaves, no_of_equilibration_sweeps, thinning_level,
        sampling_frequency)
    transposed_power_spectra = power_spectra_of_correlators[:, 1].transpose()
    norm_of_spectra_in_auxiliary_frequency_space = np.array(
        [np.absolute(item) for index, item in enumerate(np.fft.fft(transposed_power_spectra).transpose())
         if 0 < index == 2 ** (math.floor(math.log(index, 2)))])
    return [np.array([2 ** index for index in range(no_of_auxiliary_frequency_octaves + 1)]) * sampling_frequency
            / base_time_period_shift / 2 ** (no_of_auxiliary_frequency_octaves + 1), power_spectra_of_correlators[0, 0],
            norm_of_spectra_in_auxiliary_frequency_space]


def get_single_sample_of_power_trispectrum_zero_mode(algorithm_name, observable_string, output_directory,
                                                          temperature, temperature_index, no_of_sites,
                                                          no_of_auxiliary_frequency_octaves=6, base_time_period_shift=1,
                                                          no_of_equilibration_sweeps=None, thinning_level=None,
                                                          sampling_frequency=None):
    r"""
    Returns an estimate of a single sample int lim_{T -> inf}[| \Delta \tilde{Y}_T(f, s) | ** 2 / T] ds of the
    zero (auxiliary-frequency) mode S_X^3(f, f' = 0) := int lim_{T -> inf}{E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} ds
    of the power trispectrum of the time series X(t) of observable_string, where the correlator
    Y(t, s) := X(t) * X(t - s), T is the total simulation time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of
    the truncated mean-zero correlator \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is
    time, s is the auxiliary time, f is frequency, f' is the auxiliary frequency and E[.] is the expected value of the
    argument.  X(t) is considered a single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The single sample of the zero mode of the power trispectrum.  A two-dimensional numpy array of shape
        (2, N / 2 - 1) [(2, (N - 1) / 2)] for N even [odd].  Each element is a float.  The first / second sub-array is
        the frequencies / values of the zero mode of the power trispectrum.  If N is even, the frequency spectrum is
        reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum is
        reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because the correlator power
        spectra are symmetric about f = 0 and their f = 0 values are invalid for a finite-time stationary signal.
    """
    power_spectra_of_trispectrum_correlators = get_power_spectra_of_trispectrum_correlators(
        algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
        base_time_period_shift, no_of_auxiliary_frequency_octaves, no_of_equilibration_sweeps, thinning_level,
        sampling_frequency)
    return np.concatenate([np.mean(power_spectra_of_trispectrum_correlators[:, 0], axis=0),
                           np.sum(power_spectra_of_trispectrum_correlators[:, 1], axis=0)]).reshape((2, -1))


def get_single_sample_of_power_trispectrum_nonzero_mode(algorithm_name, observable_string, output_directory,
                                                             temperature, temperature_index, no_of_sites,
                                                             no_of_auxiliary_frequency_octaves=6,
                                                             base_time_period_shift=1, target_auxiliary_frequency=None,
                                                             no_of_equilibration_sweeps=None, thinning_level=None,
                                                             sampling_frequency=None):
    r"""
    Returns an estimate of a single sample
    | int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) | ds of a single
    auxiliary-frequency mode of the shortcut estimator
    E[| int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |] of the complex norm
    of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds of the
    time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total simulation
    time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, s is the auxiliary time, f is
    frequency, f' is the auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a
    single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of(\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the (correlator)
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    target_auxiliary_frequency : float or None, optional
        The target auxiliary frequency.  The auxiliary-frequency mode (of the trispectrum) whose auxiliary-frequency
        value is closest to target_auxiliary_frequency is chosen.  If None is given, the largest auxiliary frequency is
        chosen.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    List[numpy.ndarray]
        The single sample of the single mode of the power trispectrum.  A list of length 3.  The first component
        is the auxiliary frequency and is float.  The second and third components are one-dimensional numpy arrays (of
        floats) of length N / 2 - 1 [(N - 1) / 2] for N even [odd].  The second component is the frequencies.  If N is
        even, the frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is
        odd, the frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is
        because the correlator power spectra are symmetric about f = 0 and their f = 0 values are invalid for a
        finite-time stationary signal.  The third component is the trispectrum at the auxiliary-frequency value given
        by the first component.
    """
    sampling_frequency = get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_index)
    auxiliary_frequencies = (np.array([index for index in range(1, 2 ** (no_of_auxiliary_frequency_octaves + 1))])
                             * sampling_frequency / base_time_period_shift
                             / 2 ** (no_of_auxiliary_frequency_octaves + 1))
    if target_auxiliary_frequency is None:
        auxiliary_frequency_mode = 2 ** no_of_auxiliary_frequency_octaves
    else:
        auxiliary_frequency_mode = np.abs(auxiliary_frequencies - target_auxiliary_frequency).argmin()
    auxiliary_frequency = auxiliary_frequencies[auxiliary_frequency_mode]
    '''auxiliary_frequency = auxiliary_frequency_mode * sampling_frequency / base_time_period_shift / 2 ** (
            no_of_auxiliary_frequency_octaves + 1)'''
    power_spectra_of_correlators = get_power_spectra_of_trispectrum_correlators(
        algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
        base_time_period_shift, no_of_auxiliary_frequency_octaves, no_of_equilibration_sweeps, thinning_level,
        sampling_frequency)
    norm_of_spectrum_at_auxiliary_frequency_mode = np.absolute(
        np.fft.fft(power_spectra_of_correlators[:, 1].transpose()).transpose()[auxiliary_frequency_mode])
    return [auxiliary_frequency, power_spectra_of_correlators[0, 0], norm_of_spectrum_at_auxiliary_frequency_mode]


"""Base methods"""


def get_time_series(observable_string, output_directory, temperature, temperature_index, no_of_sites,
                    no_of_equilibration_sweeps=None, thinning_level=None):
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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.

    Returns
    -------
    numpy.ndarray
        The time series / sample.  A two-dimensional numpy array of shape (n, T / \Delta t), where n >= 1 (an integer)
        is the number of Cartesian components of observable_string, T is the total simulation time and \Delta t is the
        sampling interval.
    """
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    sample = get_sample_method(output_directory, temperature, temperature_index, no_of_sites,
                               no_of_equilibration_sweeps, thinning_level)
    sample = np.atleast_2d(sample)
    if len(sample) > 1:
        return sample.transpose()
    return sample


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
        (n, T / \Delta t), where n >= 1 (an integer) is the number of components of time_series and T is the total
        simulation time.
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


def get_component_averaged_power_spectrum(time_series, sampling_frequency):
    r"""
    Returns an estimate of the (Cartesian-)component average of a single sample
    lim_{T -> inf} [| \Delta \tilde{X}_T(f) | ** 2 / T] of the power spectrum
    S_X(f) := lim_{T -> inf} {E[| \Delta \tilde{X}_T(f) | ** 2] / T} of the time series X(t), where T is the total
    simulation time, \Delta \tilde{X}_T(f) is the Fourier transform of the truncated signal mean-zero time series
    \Delta X_T(t) := {X(t) - E[X] for all |t| <= T / 2, 0 otherwise}, t is time, f is frequency and E[.] is the
    expected value of the argument.  X(t) (a single component of time_series) is considered a single sample of the
    dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero time series is
    \Delta \tilde{X}_T(f_k) := sum_{n = 0}^{N − 1} \Delta X(t_n) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_X(f_k) := lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣ \Delta \tilde{X}_T(f_k) ∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series.  The factor of
    (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units of the continuum expression.

    Parameters
    ----------
    time_series : numpy.ndarray
        The time series / sample whose component-averaged single-sample power spectrum is to be computed.  A
        two-dimensional numpy array of shape (n, T / \Delta t), where n >= 1 (an integer) is the number of components
        of time_series.
    sampling_frequency : float
        The sampling frequency.

    Returns
    -------
    numpy.ndarray
        The component average of the single sample of the power spectrum.  A two-dimensional numpy array of shape
        (2, N / 2 - 1) [(2, (N - 1) / 2)] for N even [odd].  Each element is a float.  The first / second sub-array is
        the frequencies / values of the component average of the single sample of the power spectrum.  If N is
        even, the frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is
        odd, the frequency spectrum is reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is
        because the power spectrum is symmetric about f = 0 and its f = 0 values is invalid for a finite-time
        stationary signal.
    """
    power_spectra = np.atleast_2d([signal.periodogram(component_2, fs=sampling_frequency) for component_2 in
                                   [component_1 - np.mean(component_1) for component_1 in time_series]])
    # the `[:, 1:]' below removes the f = 0 value as this value is invalid for a finite-time stationary signal
    return np.mean(power_spectra, axis=0)[:, 1:]


def get_power_spectra_of_trispectrum_correlators(algorithm_name, observable_string, output_directory, temperature,
                                                 temperature_index, no_of_sites, base_time_period_shift,
                                                 no_of_auxiliary_frequency_octaves, no_of_equilibration_sweeps=None,
                                                 thinning_level=None, sampling_frequency=None):
    r"""
    Returns a numpy array of single samples of estimates of the quantity
    lim_{T -> inf}[\Delta \tilde{Y}_T(f; s) ∣ ** 2 / T] used to estimate the shortcut estimator
    E[| int lim_{T -> inf}{| \Delta \tilde{Y}_T(f; s) ∣ ** 2 / T} exp(- 2 * pi * i * f' * s) ds |] of the complex norm
    of the power trispectrum
    S_X^3(f, f') := int lim_{T -> inf} {E[| \Delta \tilde{Y}_T(f, s) | ** 2] / T} exp(- 2 * pi * i * f' * s)ds of the
    time series X(t) of observable_string, where the correlator Y(t, s) := X(t) * X(t - s), T is the total simulation
    time, \Delta \tilde{Y}_T(f, s) is the Fourier transform of the truncated mean-zero correlator
    \Delta Y_T(t, s) := {Y(t, s) - E[Y] for all |t| <= T / 2, 0 otherwise}, t is time, s is the auxiliary time, f is
    frequency, f' is the auxiliary frequency and E[.] is the expected value of the argument.  X(t) is considered a
    single sample of the dynamical 'distribution'.

    The discrete-time Fourier transform of the truncated mean-zero correlator is
    \tilde{Y}_T(f_k, s) := sum_{n = 0}^{N − 1} \Delta Y(t_n, s) exp(- 2 * pi * i * f_k * t_n), so that the estimate of
    its power spectrum is S_Y(f_k, s) := lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * (\Delta t) ** 2 / T}
    = lim_{T -> inf} {E[∣∣ \Delta \tilde{Y}_T(f_k, s) ∣∣ ** 2] * \Delta t / N}, where \Delta t is the physical sampling
    interval, t_{n + 1} = t_n + \Delta t for all n, f_k \in {0, 1 / (N \Delta t), ..., (N - 1) / (N \Delta t)} is the
    discrete frequency spectrum and N = T / \Delta t is the number of samples of the time series (of the
    correlator).  The factor of (\Delta t) ** 2 in the definition of the discrete-time power spectrum retains the units
    of the continuum expression.

    If observable_string is a Cartesian vector of dimension greater than 1, each single sample of the quantity is
    computed for each Cartesian component and the average of the resultant quantities are returned.

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
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_auxiliary_frequency_octaves : int, optional
        The number of auxiliary-frequency octaves over which the trispectrum is estimated.
    base_time_period_shift : int, optional
        The elementary number of multiples of the physical sampling interval by which the time series is shifted in
        order to form the correlator whose power spectrum is computed in order to compute the trispectrum.
    no_of_equilibration_sweeps : None or int, optional
        The number of discarded equilibration samples.  If None, no_of_equilibration_sweeps = 0.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.
    sampling_frequency : float or None, optional
        The sampling frequency.  If None is given, a float is computed via sample_getter.get_physical_time_step(),
        which computes the emergent physical time step of the Metropolis algorithm, where the physical timescale arises
        due to the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Returns
    -------
    numpy.ndarray
        The (single samples of the) power spectra of the trispectrum correlators.  A two-dimensional numpy array
        of shape (no_of_auxiliary_frequency_octaves + 1, m / 2 - 1)
        [(no_of_auxiliary_frequency_octaves + 1, (m - 1) / 2)] for m  even [odd], where
        m = N - (2 ** (no_of_auxiliary_frequency_octaves + 1) - 1) * base_time_period_shift.  The first / second
        sub-array is the frequencies / values of the (single samples of the) power spectra.  If m is even, the
        frequency spectrum is reduced to f_k \in {1 / (m \Delta t), ..., (m / 2 - 1) / (m \Delta t)}; if m is odd, the
        frequency spectrum is reduced to f_k \in {1 / (m \Delta t), ..., (m - 1) / 2 / (m \Delta t)}.  This is because
        the correlator power spectra are symmetric about f = 0 and their f = 0 values are invalid for a finite-time
        stationary signal.
    """
    sampling_frequency = get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_index)
    time_series = get_time_series(observable_string, output_directory, temperature, temperature_index, no_of_sites,
                                  no_of_equilibration_sweeps, thinning_level)
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


def get_power_spectrum_means_and_std_errors(power_spectra):
    r"""
    Returns both 1) the mean of the power spectra estimated from M repeated simulations, and 2) an estimate the
    standard error (of the spectrum) with respect to the repeated simulations.

    Parameters
    ----------
    power_spectra : List[numpy.ndarray]
        The power spectra estimated from M repeated simulations.  A list of length M.  Component i is the power
        spectrum estimated from simulation i and is a numpy array (of floats) of length N / 2 - 1 [(N - 1) / 2] for N
        even [odd], where N is the number of samples within each simulation.  If N is even, the frequency spectrum
        is reduced to f_k \in {1 / (N \Delta t), ..., (N / 2 - 1) / (N \Delta t)}; if N is odd, the frequency spectrum
        is reduced to f_k \in {1 / (N \Delta t), ..., (N - 1) / 2 / (N \Delta t)}.  This is because each power spectrum
        is symmetric about f = 0 and its f = 0 value is invalid for a finite-time stationary signal.

    Returns
    -------
    numpy.ndarray
        The power spectrum frequencies, means and standard errors.  A two-dimensional numpy array of shape (3, N).  The
        standard error at frequency f is the square root of the ratio of the simulation variance to M, where the sample
        is the dynamical sample (at f) whose samples are the repeated simulations.
    """
    power_spectra = np.array(power_spectra)
    return np.concatenate([np.mean(power_spectra, axis=0).flatten(),
                           np.std(power_spectra, axis=0)[1] / len(power_spectra) ** 0.5]).reshape((3, -1))
