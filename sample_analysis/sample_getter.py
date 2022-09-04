"""This module contains methods that convert the entire sample [stored (by the Fortran code) in the sample.csv file at
    each temperature and repeated job] into sample observables.  The number of ECMC events and rate of Metropolis
    acceptances are retrieved by other methods; the sampling frequency and emergent-Langevin physical timescale are
    computed by further methods."""
import math
import numpy as np


"""Methods common to all models"""


def get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_index):
    """
    Returns the physical sampling frequency of the Metropolis algorithms, where the physical timescale arises due to
    the diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Parameters
    ----------
    algorithm_name : str
        The name of the algorithm used to generate the sample.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    sampling_frequency : float or None
        The emergent physical sampling frequency.  If None is given, a float is computed via get_physical_time_step();
        if a float is given, the float is returned.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    float
        The emergent physical sampling frequency.
    """
    if sampling_frequency is None:
        return 1.0 / get_physical_time_step(algorithm_name, output_directory, temperature_index)
    return sampling_frequency


def get_physical_time_step(algorithm_name, output_directory, temperature_index):
    """
    Returns the physical time step of the Metropolis algorithms, where the physical timescale arises due to the
    diffusive Langevin dynamics that emerges from Metropolis dynamics.

    Parameters
    ----------
    algorithm_name : str
        The name of the algorithm used to generate the sample.
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    float
        The emergent physical time step.
    """
    acceptance_rates = get_acceptance_rates(output_directory, temperature_index)
    if algorithm_name == "hxy-gaussian-noise-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis":
        return 0.5 * acceptance_rates[1] * acceptance_rates[0] ** 2  # emergent Brownian diffusivity for Gaussian noise
    return acceptance_rates[1] * acceptance_rates[0] ** 2 / 24.0  # emergent Brownian diffusivity for uniform noise


def get_acceptance_rates(output_directory, temperature_index):
    """
    Returns the acceptance rates of the Metropolis algorithms.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    numpy.ndarray
        The acceptance rates.  A one-dimensional numpy array of length 2, 3 or 4.  Each element is a float.  The
        first / second element is the final width of the proposal interval for the rotational moves ('final' since this
        width can be tuned during equilibration in order to reach some target acceptance rate) / the acceptance rate of
        the field rotations.  If the length of the array is 3 and the model is an XY model, the third element is the
        acceptance rate of the external global moves (the externally applied global spin twists).  If the model is a
        Maggs-electrolyte model, the third element is the acceptance rate of the charge hops.  If the length of the
        array is 4, the model is a Maggs-electrolyte model and the fourth element is the acceptance rate of the
        external global moves (the externally applied topological-sector fluctuations).
    """
    return np.atleast_1d(np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/acceptance_rates.csv",
                                    dtype=float, delimiter=","))


def get_no_of_events(output_directory, temperature_index):
    """
    Returns the number of ECMC events.  If external global Metropolis moves (the externally applied global spin twists)
    are chosen, this method also returns the acceptance rates of these moves.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    numpy.ndarray
        The number of events and the acceptance rate of the external global Metropolis moves (if chosen).  A
        one-dimensional numpy array of length 1 or 2.  The first element is the number of ECMC events and has type int.
        If the length of the array is 2, external global Metropolis moves (the global spin twists) have been chosen and
        the second element is the acceptance rate of these moves and has type float.
    """
    return np.atleast_1d(np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/no_of_events.csv", dtype=float,
                                    delimiter=","))


def get_potential(output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps=None,
                  thinning_level=None):
    """
    Returns the potential sample.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The potential sample.  A one-dimensional numpy array of length no_of_observations.  The nth element is a float
        corresponding to the potential measured at observation n.
    """
    return get_reduced_sample(np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/potential.csv", dtype=float,
                                         delimiter=","), no_of_equilibration_sweeps, thinning_level)


def get_specific_heat(output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps=None,
                      thinning_level=None):
    """
    Returns the sample of the specific-heat C_V(x; temperature, no_of_sites) per particle, where
    C_V(x; temperature, no_of_sites) = [U(x; temperature, no_of_sites) - E[U(x; temperature, no_of_sites)]] ** 2 /
    temperature ** 2, with x the particle positions at the time of observation, E[.] the expected value of the argument
    and U(x; temperature, no_of_sites) the system potential.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The specific-heat sample.  A one-dimensional numpy array of length no_of_observations.  The nth element is a
        float corresponding to the specific heat measured at observation n.
    """
    potential_sample = get_potential(output_directory, temperature, temperature_index, no_of_sites,
                                     no_of_equilibration_sweeps, thinning_level)
    return (potential_sample - np.mean(potential_sample)) ** 2 / no_of_sites / temperature ** 2


def get_external_global_move(output_directory, temperature, temperature_index, no_of_sites,
                             no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the external global move vector, which is an integer-valued two-dimensional Cartesian vector
        whose x / y component is 0 if the x / y component of the external global move was not accepted, 1 if a positive
        external global move was accepted and -1 if a negative external global move was accepted.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The external global move sample.  A two-dimensional numpy array of shape (no_of_observations, 2).  The nth
        sub-array is an integer-valued two-dimensional Cartesian vector, as measured at observation n; the first /
        second element (of the sub-array) is an int corresponding to the x / y Cartesian component, where the value is
        0 if the x / y component of the external global move was not accepted, 1 if a positive external global move was
        accepted and -1 if a negative external global move was accepted.
    """
    return get_reduced_sample(np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/external_global_moves.csv",
                                         dtype=int, delimiter=","), no_of_equilibration_sweeps, thinning_level)


"""XY and HXY magnetisation methods"""


def get_non_normalised_cartesian_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                               no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the non-normalised magnetisation no_of_sites * m(x; temperature, no_of_sites), where
    m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian magnetisation, with
    x_i the position of particle i at the time of observation.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the non-normalised magnetisation.  A two-dimensional numpy array of shape (no_of_observations, 2).
        The nth sub-array is the non-normalised magnetisation vector measured at observation n; its first / second
        element is a float corresponding to the x / y component of the non-normalised magnetisation vector measured at
        observation n.
    """
    return get_reduced_sample(
        np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/magnetisation.csv", dtype=float, delimiter=","),
        no_of_equilibration_sweeps, thinning_level)


def get_magnetisation_norm(output_directory, temperature, temperature_index, no_of_sites,
                           no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the magnetisation norm || m(x; temperature, no_of_sites) ||, where
    m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian magnetisation, with
    x_i the position of particle i at the time of observation.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the magnetisation norm.  A one-dimensional numpy array of length no_of_observations.  The nth
        element is a float corresponding to magnetisation norm measured at observation n.
    """
    return np.linalg.norm(get_non_normalised_cartesian_magnetisation(
        output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps, thinning_level),
        axis=1) / no_of_sites


def get_magnetisation_phase(output_directory, temperature, temperature_index, no_of_sites,
                            no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the magnetisation phase phi(x; temperature, no_of_sites), where
    m(x; temperature, no_of_sites) = (|| m(x; temperature, no_of_sites) ||, phi(x; temperature, no_of_sites))^t in
    radial coordinates and m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian
    magnetisation, with x_i the position of particle i at the time of observation.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the magnetisation phase.  A one-dimensional numpy array of length no_of_observations.  The nth
        element is a float corresponding to magnetisation phase measured at observation n.
    """
    return np.array([get_phase_in_polar_coordinates(observation) for observation in
                     get_non_normalised_cartesian_magnetisation(output_directory, temperature, temperature_index,
                                                                no_of_sites, no_of_equilibration_sweeps,
                                                                thinning_level)])


def get_rotated_magnetisation_phase(output_directory, temperature, temperature_index, no_of_sites,
                                    no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the rotated magnetisation phase, where the magnetisation phase
    phi(x; temperature, no_of_sites) is defined by writing the magnetisation vector
    m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites in radial coordinates,
    m(x; temperature, no_of_sites) = (|| m(x; temperature, no_of_sites) ||, phi(x; temperature, no_of_sites))^t.  Here,
    x_i the position of particle i at the time of observation.

    Each observation of the phase is rotated by the mean of the absolute values of the observations.  We then apply
    modular arithmetic to transform the resultant values onto (-pi, pi].  We then subtract the mean of the resultant
    sample.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the rotated magnetisation phase.  A one-dimensional numpy array of length no_of_observations.
        The nth element is a float corresponding to the rotated magnetisation phase measured at observation n.
    """
    non_rotated_mag_phase = get_magnetisation_phase(output_directory, temperature, temperature_index, no_of_sites,
                                                    no_of_equilibration_sweeps, thinning_level)
    return (non_rotated_mag_phase - np.sign(np.mean(non_rotated_mag_phase)) * np.mean(abs(non_rotated_mag_phase))
            + math.pi) % (2.0 * math.pi) - math.pi


def get_magnetisation_squared(output_directory, temperature, temperature_index, no_of_sites,
                              no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the magnetisation squared || m(x; temperature, no_of_sites) ||^2, where
    m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian magnetisation, with
    x_i the position of particle i at the time of observation.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the magnetisation norm.  A one-dimensional numpy array of length no_of_observations.  The nth
        element is a float corresponding to magnetisation norm measured at observation n.
    """
    return get_magnetisation_norm(output_directory, temperature, temperature_index, no_of_sites,
                                  no_of_equilibration_sweeps, thinning_level) ** 2


def get_cartesian_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the Cartesian magnetisation vector m(x; temperature, no_of_sites) =
    sum_i [cos(x_i), sin(x_i)]^t / no_of_sites, where x_i the position of particle i at the time of observation.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the Cartesian magnetisation.  A two-dimensional numpy array of shape (no_of_observations, 2).
        The nth sub-array is the Cartesian magnetisation vector measured at observation n; its first / second element
        is a float corresponding to the x / y component of the Cartesian magnetisation vector measured at observation n.
    """
    return get_non_normalised_cartesian_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                                      no_of_equilibration_sweeps, thinning_level) / no_of_sites


def get_absolute_cartesian_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                         no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the vector formed from the absolute value of each Cartesian component of the magnetisation
    vector m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites, where x_i the position of
    particle i at the time of observation.  Each observation is therefore the vector
    [|m_1(x; temperature, no_of_sites)|, |m_2(x; temperature, no_of_sites)|].

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the absolute Cartesian magnetisation vector.  A two-dimensional numpy array of shape
        (no_of_observations, 2).  The nth sub-array is the absolute Cartesian magnetisation vector measured at
        observation n; its first / second element is a float corresponding to the absolute value of the x / y component
        of the Cartesian magnetisation vector measured at observation n.
    """
    return np.abs(get_cartesian_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                              no_of_equilibration_sweeps, thinning_level))


def get_magnetic_susceptibility(output_directory, temperature, temperature_index, no_of_sites,
                                no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the magnetic-norm susceptibility per particle chi_m(x; temperature, no_of_sites) =
    no_of_sites * [|| m(x; temperature, no_of_sites) || - E[|| m(x; temperature, no_of_sites) ||]] ** 2 / temperature,
    where m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian magnetisation,
    with x_i the position of particle i at the time of observation.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the magnetic susceptibility.  A one-dimensional numpy array of length no_of_observations.  The nth
        element is a float corresponding to the magnetisation phase measured at observation n.
    """
    magnetisation_norm = get_magnetisation_norm(output_directory, temperature, temperature_index, no_of_sites,
                                                no_of_equilibration_sweeps, thinning_level)
    return no_of_sites * (magnetisation_norm - np.mean(magnetisation_norm)) ** 2 / temperature


def get_cartesian_relative_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                         no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the Cartesian relative magnetisation vector tilde{m}(x; temperature, no_of_sites) =
    m / sigma_{|| m ||}, where m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the
    Cartesian magnetisation (with x_i the position of particle i at the time of observation) and sigma_{|| m ||} is
    the standard deviation of its norm.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the Cartesian relative magnetisation.  A two-dimensional numpy array of shape
        (no_of_observations, 2).  The nth sub-array is the Cartesian magnetisation vector measured at observation n;
        its first / second element is a float corresponding to the x / y component of the relative Cartesian
        magnetisation vector measured at observation n.
    """
    cartesian_magnetisation = get_cartesian_magnetisation(output_directory, temperature, temperature_index, no_of_sites,
                                                          no_of_equilibration_sweeps, thinning_level)
    return cartesian_magnetisation / np.std(np.linalg.norm(cartesian_magnetisation, axis=1))


def get_relative_magnetisation_norm(output_directory, temperature, temperature_index, no_of_sites,
                                    no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the relative magnetisation norm || tilde{m}(x; temperature, no_of_sites) ||, where
    tilde{m}(x; temperature, no_of_sites) = m / sigma_{|| m ||} is the Cartesian relative magnetisation,
    m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian magnetisation (with
    x_i the position of particle i at the time of observation) and sigma_{|| m ||} is the standard deviation of its
    norm.


    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the relative magnetisation norm.  A one-dimensional numpy array of length no_of_observations.
        The nth element is a float corresponding to the relative magnetisation norm measured at observation n.
    """
    return np.linalg.norm(get_cartesian_relative_magnetisation(output_directory, temperature, temperature_index,
                                                               no_of_sites, no_of_equilibration_sweeps,
                                                               thinning_level), axis=1)


"""XY and HXY emergent-field-based methods"""


def get_non_normalised_total_vortex_polarisation(output_directory, temperature, temperature_index, no_of_sites,
                                                 no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the non-normalised total vortex polarisation vector, where the x / y component of the total
    vortex polarisation is the normalised (with respect to no_of_sites) sum over the first derivatives of the potential
    with respect to the spin differences along the x / y dimension.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the non-normalised total vortex polarisation.  A two-dimensional numpy array of shape
        (no_of_observations, 2).  The nth sub-array is the non-normalised total vortex polarisation measured at
        observation n; its first / second element is a float corresponding to the x / y component of the non-normalised
        total vortex polarisation measured at observation n.
    """
    return get_reduced_sample(np.loadtxt(
        f"{output_directory}/temp_{temperature_index:02d}/1st_deriv_of_potential.csv", dtype=float, delimiter=","),
        no_of_equilibration_sweeps, thinning_level)


def get_inverse_vacuum_permittivity(output_directory, temperature, temperature_index, no_of_sites,
                                    no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the inverse vacuum permittivity vector, whose x / y component is the normalised (with respect
    to no_of_sites) sum over the second derivatives of the potential with respect to the spin differences along the
    x / y dimension.  This is a generalisation of its HXY definition in J. Phys. Condens. Matter 29, 085402 (2017).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the inverse vacuum permittivity.  A one-dimensional numpy array of length no_of_observations.
        The nth  element is a float corresponding to the inverse vacuum permittivity measured at observation n.
    """
    return get_reduced_sample(np.mean(np.loadtxt(
        f"{output_directory}/temp_{temperature_index:02d}/2nd_deriv_of_potential.csv", dtype=float, delimiter=","),
        axis=1), no_of_equilibration_sweeps, thinning_level) / no_of_sites


def get_total_vortex_polarisation(output_directory, temperature, temperature_index, no_of_sites,
                                  no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the total vortex polarisation vector, whose x / y component is the normalised (with respect
    to no_of_sites) sum over the first derivatives of the potential with respect to the spin differences along the
    x / y dimension.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the total vortex polarisation.  A two-dimensional numpy array of shape (no_of_observations, 2).
        The nth sub-array is the total vortex polarisation measured at observation n; its first / second element is a
        float corresponding to the x / y component of the total vortex polarisation measured at observation n.
    """
    return get_non_normalised_total_vortex_polarisation(output_directory, temperature, temperature_index,
                                                        no_of_sites, no_of_equilibration_sweeps,
                                                        thinning_level) / no_of_sites


def get_helicity_modulus(output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps=None,
                         thinning_level=None):
    r"""
    Returns the sample of the helicity modulus Upsilon(x; temperature, no_of_sites).  For the XY model,
    Upsilon(x; temperature, no_of_sites) = sum_{<i, j>} cos(x_i - x_j) / no_of_sites - [sum_{<i, j>} sin(x_i - x_j)
    - E[sum_{<i, j>} sin(x_i - x_j)]] ** 2 / temperature / no_of_sites.  For the HXY model,
    Upsilon(x; temperature, no_of_sites) = sum_{n = 1}^{inf} sum_{<i, j>} (-1) ** n * cos[n * (x_i - x_j)]
    / no_of_sites - no_of_sites * [\overline{E} - E[\overline{E}]] ** 2 / temperature.  In both equations, E[.] is the
    expected value of the argument; in the latter, \overline{E} is the zero mode of the emergent electric field, as
    defined in J. Phys. Condens. Matter 29, 085402 (2017).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the helicity modulus.  A one-dimensional numpy array of length no_of_observations.  The nth
        element is a float corresponding to the helicity modulus measured at observation n.
    """
    non_normalised_total_vortex_polarisation = get_non_normalised_total_vortex_polarisation(
        output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps, thinning_level)
    return (get_inverse_vacuum_permittivity(output_directory, temperature, temperature_index, no_of_sites,
                                            no_of_equilibration_sweeps, thinning_level)
            - np.mean((non_normalised_total_vortex_polarisation
                       - np.mean(non_normalised_total_vortex_polarisation, axis=0)) ** 2, axis=1)
            / temperature / no_of_sites)


def get_potential_minimising_twists(output_directory, temperature, temperature_index, no_of_sites,
                                    no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the potential-minimising twist field -- an integer-valued two-dimensional vector field whose
    x / y component corresponds to the number of twists required (along the x / y dimension) to minimise the potential.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the potential-minimising twist field.  A two-dimensional numpy array of shape
        (no_of_observations, 2).  The nth sub-array is the potential-minimising twist field measured at observation n;
        its first / second element is an int corresponding to the x / y component of the potential-minimising twist
        field measured at observation n.
    """
    return get_reduced_sample(np.loadtxt(
        f"{output_directory}/temp_{temperature_index:02d}/potential_minimising_twists.csv", dtype=float, delimiter=","),
        no_of_equilibration_sweeps, thinning_level)


def get_potential_minimising_twist_susceptibility(output_directory, temperature, temperature_index, no_of_sites,
                                                  no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the potential-minimising twist susceptibility chi_{\tilde{t}}(x; temperature, no_of_sites) =
    4 * \pi ** 2 * [\tilde{t} - E[\tilde{t}]] ** 2 / temperature, where E[.] is the expected value of the
    argument and \tilde{t} \in Z^2 is the potential-minimising twist field w -- an integer-valued two-dimensional
    vector field whose x / y component corresponds to the number of twists required (along the x / y dimension) to
    minimise the potential.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the potential-minimising twist susceptibility.  A one-dimensional numpy array of length
        no_of_observations.  The nth element is a float corresponding to the potential-minimising twist susceptibility
        measured at observation n.
    """
    potential_minimising_twists_sample = get_potential_minimising_twists(
        output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps, thinning_level)
    return 4.0 * math.pi ** 2 * np.mean(
        (potential_minimising_twists_sample - np.mean(potential_minimising_twists_sample, axis=0)) ** 2,
        axis=1) / temperature


def get_hxy_topological_sector(output_directory, temperature, temperature_index, no_of_sites,
                               no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the HXY topological sector w \in Z^2, an integer-valued vector field whose components are
    fixed such that the minimal toroidal vortex polarisation vector \overline{E}_{p, x / y} \in ( - \pi / L, \pi / L],
    where \overline{E} = \overline{E}_p + 2 \pi w / L is the zero mode of the emergent electric field, as defined in
    J. Phys. Condens. Matter 29, 085402 (2017).  This method is only valid for the HXY model.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the topological sector.  A two-dimensional numpy array of shape (no_of_observations, 2).  The nth
        sub-array is the topological sector measured at observation n; its first / second element is an int
        corresponding to the x / y component of the topological sector measured at observation n.
    """
    non_normalised_total_vortex_polarisation = get_non_normalised_total_vortex_polarisation(
        output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps, thinning_level)
    return (non_normalised_total_vortex_polarisation + math.pi * no_of_sites ** 0.5) // (2.0 * math.pi *
                                                                                         no_of_sites ** 0.5)


def get_hxy_topological_susceptibility(output_directory, temperature, temperature_index, no_of_sites,
                                       no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the HXY topological susceptibility chi_w^H(x; temperature, no_of_sites) =
    no_of_sites * [\overline{E}_w - E[\overline{E}_w]] ** 2 / temperature, where E[.] is the expected value of the
    argument and \overline{E}_w = 2 \pi w / L is the topological sector w \in Z^2, an integer-valued vector field whose
    components are fixed such that the minimal toroidal vortex polarisation vector
    \overline{E}_{p, x / y} \in ( - \pi / L, \pi / L], where \overline{E} = \overline{E}_p + 2 \pi w / L is the zero
    mode of the emergent electric field, as defined in J. Phys. Condens. Matter 29, 085402 (2017).  This method is only
    valid for the HXY model.

    In Phys. Rev. B 91, 155412 (2015), the HXY topological susceptibility was called the HXY winding-field
    susceptibility and was defined without the factor of 1/2 to account for each Cartesian dimension of the system; in
    the return line below, this corresponds to the first np.mean() -> np.sum().

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the topological susceptibility.  A one-dimensional numpy array of length no_of_observations.  The
        nth element is a float corresponding to the topological susceptibility measured at observation n.
    """
    topological_sector_sample = get_hxy_topological_sector(output_directory, temperature, temperature_index,
                                                           no_of_sites, no_of_equilibration_sweeps, thinning_level)
    return 4.0 * math.pi ** 2 * np.mean((topological_sector_sample - np.mean(topological_sector_sample, axis=0)) ** 2,
                                        axis=1) / temperature


"""Maggs-electrolyte methods"""


def get_sum_of_electric_field(output_directory, temperature, temperature_index, no_of_sites,
                              no_of_equilibration_sweeps=None, thinning_level=None):
    """
    Returns the sample of the sum (over the lattice sites) of the electric field.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the sum of the electric field.  A two-dimensional numpy array of shape (no_of_observations, 2).
        The nth sub-array is the sum of the electric field measured at observation n; its first / second element is a
        float corresponding to the x / y component of the sum of the electric field measured at observation n.
    """
    return get_reduced_sample(np.loadtxt(
        f"{output_directory}/temp_{temperature_index:02d}/electric_field_sum.csv", dtype=float, delimiter=","),
        no_of_equilibration_sweeps, thinning_level)


def get_electric_field_zero_mode(output_directory, temperature, temperature_index, no_of_sites,
                                 no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the zero mode of the electric field \overline{E} = sum_i E(i) / no_of_sites, where E(i) is
    the electric field at lattice site i.  The zero mode was defined as the `harmonic mode' in
    Phys. Rev. B 91, 155412 (2015).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the zero mode of the electric field.  A two-dimensional numpy array of shape
        (no_of_observations, 2).  The nth sub-array is the zero mode of the electric field measured at observation n;
        its first / second element is a float corresponding to the x / y component of the zero mode of the electric
        field measured at observation n.
    """
    return get_sum_of_electric_field(output_directory, temperature, temperature_index, no_of_sites,
                                     no_of_equilibration_sweeps, thinning_level) / no_of_sites


def get_inverse_permittivity(output_directory, temperature, temperature_index, no_of_sites,
                             no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the inverse (electric) permittivity modulus [epsilon(x; temperature, no_of_sites)] ** (-1)
    = 1.0 - no_of_sites * [\overline{E} - E[\overline{E}]] ** 2 / temperature, where E[.] is the expected value of the
    argument and \overline{E} is the zero mode of the electric field, as defined in Phys. Rev. B 91, 155412 (2015).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the helicity modulus.  A one-dimensional numpy array of length no_of_observations.  The nth
        element is a float corresponding to the helicity modulus measured at observation n.
    """
    sum_of_electric_field_sample = get_sum_of_electric_field(output_directory, temperature, temperature_index,
                                                             no_of_sites, no_of_equilibration_sweeps, thinning_level)
    return 1.0 - np.mean((sum_of_electric_field_sample - np.mean(sum_of_electric_field_sample, axis=0)) ** 2,
                         axis=1) / temperature / no_of_sites


def get_topological_sector(output_directory, temperature, temperature_index, no_of_sites,
                           no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the topological sector w \in Z^2, an integer-valued vector field whose components are fixed
    such that the minimal toroidal polarisation vector \overline{E}_{p, x / y} \in ( - \pi / L, \pi / L], where
    \overline{E} = \overline{E}_p + 2 \pi w / L is the zero mode of the electric field, as defined in
    Phys. Rev. B 91, 155412 (2015).  This is only strictly defined for a charge-neutral electrolyte of particles with
    charge \pm 2.0 * pi.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the topological sector.  A two-dimensional numpy array of shape (no_of_observations, 2).  The nth
        sub-array is the topological sector measured at observation n; its first / second element is an int
        corresponding to the x / y component of the topological sector measured at observation n.
    """
    sum_of_electric_field = get_sum_of_electric_field(output_directory, temperature, temperature_index, no_of_sites,
                                                      no_of_equilibration_sweeps, thinning_level)
    return (sum_of_electric_field + math.pi * no_of_sites ** 0.5) // (2.0 * math.pi * no_of_sites ** 0.5)


def get_topological_susceptibility(output_directory, temperature, temperature_index, no_of_sites,
                                   no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns the sample of the topological susceptibility chi_w(x; temperature, no_of_sites) =
    no_of_sites * [\overline{E}_w - E[\overline{E}_w]] ** 2 / temperature, where E[.] is the expected value of the
    argument and \overline{E}_w = 2 \pi w / L is the winding field, as defined in Phys. Rev. B 91, 155412 (2015).  In
    Phys. Rev. B 91, 155412 (2015), topological susceptibility was called the winding-field susceptibility and was
    defined without the factor of 1/2 to account for each Cartesian dimension of the system; in the return line below,
    this corresponds to the first np.mean() -> np.sum().

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The total number of lattice sites.
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of observations to be discarded between retained observations of the thinning process.  If None,
        all observations are retained.

    Returns
    -------
    numpy.ndarray
        The sample of the topological susceptibility.  A one-dimensional numpy array of length no_of_observations.  The
        nth element is a float corresponding to the topological susceptibility measured at observation n.
    """
    topological_sector_sample = get_topological_sector(output_directory, temperature, temperature_index, no_of_sites,
                                                       no_of_equilibration_sweeps, thinning_level)
    return 4.0 * math.pi ** 2 * np.mean((topological_sector_sample - np.mean(topological_sector_sample, axis=0)) ** 2,
                                        axis=1) / temperature


"""helper methods"""


def get_reduced_sample(entire_sample, no_of_equilibration_sweeps=None, thinning_level=None):
    if no_of_equilibration_sweeps is None:
        if thinning_level is None:
            return entire_sample
        return entire_sample[::thinning_level]
    if thinning_level is None:
        return entire_sample[no_of_equilibration_sweeps:]
    return entire_sample[no_of_equilibration_sweeps::thinning_level]


def get_phase_in_polar_coordinates(two_dimensional_cartesian_vector):
    """
    Returns the phase phi of the transformation of two_dimensional_cartesian_vector into polar coordinates.

    Parameters
    ----------
    two_dimensional_cartesian_vector : numpy.ndarray
        A two-dimensional vector in Cartesian coordinates.  A one-dimensional numpy array of length 2.  The 1st / 2nd
        element is a float corresponding to the x / y component of the Cartesian vector.

    Returns
    -------
    float
        The phase phi of the transformation of two_dimensional_cartesian_vector into polar coordinates.
    """
    if two_dimensional_cartesian_vector[0] > 0.0:
        return math.atan(two_dimensional_cartesian_vector[1] / two_dimensional_cartesian_vector[0])
    elif two_dimensional_cartesian_vector[1] > 0.0:
        return math.atan(two_dimensional_cartesian_vector[1] / two_dimensional_cartesian_vector[0]) + math.pi
    else:
        return math.atan(two_dimensional_cartesian_vector[1] / two_dimensional_cartesian_vector[0]) - math.pi
