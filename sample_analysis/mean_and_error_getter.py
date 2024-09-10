"""This module contains methods that convert the summary statistics [stored (by the Fortran code) at each temperature
    and repeated run] into the summary statistics of observables of interest."""
import numpy as np


"""Methods common to all models"""


def get_potential_summary_stats(output_directory, temperature_index):
    """
    Returns the summary stats of the potential.

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
        The summary stats of the potential.  A one-dimensional numpy array of floats of length 4.
    """
    return np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/potential_summary_stats.csv",
                      dtype=float, delimiter=",")


def get_potential(output_directory, temperature_index):
    """
    Returns the Monte Carlo mean and error of the potential.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_potential_summary_stats(output_directory, temperature_index)
    return summary_stats[0], summary_stats[1]


def get_specific_heat(output_directory, temperature_index):
    """
    Returns the Monte Carlo mean and error of the specific heat per particle
        C_V(x; temperature, no_of_sites) / no_of_sites = [U(x; temperature, no_of_sites) -
                                                        \bar{U(x; temperature, no_of_sites)}] ** 2 / temperature ** 2,
        with x the particle positions at the time of sample and U(x; temperature, no_of_sites) the system potential.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_potential_summary_stats(output_directory, temperature_index)
    return summary_stats[2], summary_stats[3]


"""XY and HXY magnetisation methods"""


def get_magnetic_summary_stats(output_directory, temperature_index):
    """
    Returns the summary stats related to the magnetisation.

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
        The summary stats of the magnetisation.  A one-dimensional numpy array of floats of length 4.
    """
    return np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/magnetic_summary_stats.csv",
                      dtype=float, delimiter=",")


def get_magnetisation_norm(output_directory, temperature_index):
    """
    Returns the Monte Carlo mean and error of the norm of the magnetisation vector.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_magnetic_summary_stats(output_directory, temperature_index)
    return summary_stats[0], summary_stats[1]


def get_magnetic_susceptibility(output_directory, temperature_index):
    """
    Returns the Monte Carlo mean and error of the magnetic-norm susceptibility per particle
        chi_m(x; temperature, no_of_sites) / no_of_sites =
                [|| m(x; temperature, no_of_sites) || - \bar{|| m(x; temperature, no_of_sites) ||}] ** 2 / temperature,
    where m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian magnetisation,
    with x_i the position of particle i at the time of sample.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
            The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st
            element is the Monte Carlo mean / error.
    """
    summary_stats = get_magnetic_summary_stats(output_directory, temperature_index)
    return summary_stats[2], summary_stats[3]


def get_helicity_modulus(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the helicity modulus.  For the XY model,
    Upsilon(x; temperature, no_of_sites) = 0.5 * sum_{<i, j>} cos(x_i - x_j) / no_of_sites - 0.5 * ([sum_{<i, j>}_x
    sin(x_i - x_j) - E[sum_{<i, j>}_x sin(x_i - x_j)]] ** 2 + [sum_{<i, j>}_y sin(x_i - x_j) -
    E[sum_{<i, j>}_y sin(x_i - x_j)]] ** 2) / temperature / no_of_sites.  For the HXY model,
    Upsilon(x; temperature, no_of_sites) = 0.5 * sum_{n = 1}^{inf} sum_{<i, j>} (-1) ** n * cos[n * (x_i - x_j)]
    / no_of_sites - 0.5 * no_of_sites * [\overline{E} - E[\overline{E}]] ** 2 / temperature.  In both equations, E[.]
    is the expected value of the argument; in the latter, \overline{E} is the zero mode of the emergent electric field,
    as defined in J. Phys.: Condens. Matter 29, 085402 (2017).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/helicity_summary_stats.csv",
                               dtype=float, delimiter=",")
    return summary_stats[0], summary_stats[1]


def get_hot_twist_relaxation_susceptibility(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the hot global twist-relaxation susceptibility
    chi_{\tilde{t}}(x; temperature, no_of_sites) = 0.5 * 4 * \pi ** 2 * [\tilde{t} - E[\tilde{t}]] ** 2 / temperature,
    where E[.] is the expected value of the argument and \tilde{t} \in Z^2 is the hot global twist-relaxation field --
    an integer-valued two-dimensional vector field whose x / y component corresponds to the number of twists required
    (along the x / y dimension) to minimise the non-annealed (ie, hot) potential.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = np.loadtxt(
        f"{output_directory}/temp_{temperature_index:02d}/hot_twist_relaxations_summary_stats.csv", dtype=float,
        delimiter=",")
    return summary_stats[0], summary_stats[1]


def get_emergent_field_summary_stats(output_directory, temperature_index):
    """
    Returns the summary stats related to the zero mode of the emergent electric field.

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
        The summary stats of the emergent electric field.  A one-dimensional numpy array of floats of length 4.
    """
    return np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/emergent_field_summary_stats.csv",
                      dtype=float, delimiter=",")


def get_emergent_field_zero_mode_susceptibility(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the XY/HXY zero-mode susceptibility
    chi_\overline{E}(x; temperature, no_of_sites) =
    0.5 * no_of_sites * [\overline{E} - E[\overline{E}]] ** 2 / temperature, where E[.] is the expected value of the
    argument.  As described in sample_getter.get_xy_emergent_field_zero_mode(), the XY zero mode is a fictitious
    quantity in some sense.

    In J. Phys.: Condens. Matter 29, 085402 (2017), the HXY zero-mode susceptibility was called the HXY harmonic-mode
    susceptibility and was defined without the factor of 1/2 to account for each Cartesian dimension of the system.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_emergent_field_summary_stats(output_directory, temperature_index)
    return summary_stats[0], summary_stats[1]


def get_global_defect_susceptibility(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the XY/HXY global-defect susceptibility
    chi_w(x; temperature, no_of_sites) = 0.5 * no_of_sites * [\overline{E}_w - E[\overline{E}_w]] ** 2 / temperature,
    where E[.] is the expected value of the argument and \overline{E}_w = 2 \pi w / L is the XY/HXY topological sector
    w \in Z^2, an integer-valued vector field defined in get_xy_topological_sector().  As described in
    sample_getter.get_xy_topological_sector(), the XY global topological defects are fictitious quantities in some
    sense.

    In the HXY case, the global-defect field is equivalent to the internal-twist field, so the HXY global-defect
    susceptibility is also called the HXY internal-twist susceptibility.

    Note that this Python method is only valid for the two-dimensional XY models.

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_emergent_field_summary_stats(output_directory, temperature_index)
    return summary_stats[2], summary_stats[3]


def get_xy_twist_relaxation_susceptibility(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the XY global twist-relaxation susceptibility, where the XY global
    twist-relaxation field is a two-dimensional integer-valued vector field whose x / y component corresponds to the
    number of externally applied global spin twists required (along the x / y dimension) to minimise the XY-annealed
    potential (ie, once spin waves (of the configuration in question) have been annealed away, keeping the topological
    defects fixed in position -- as described in Phys. Rev. B 109, 085405 (2024)).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/twist_relaxations_summary_stats.csv",
                               dtype=float, delimiter=",")
    return summary_stats[0], summary_stats[1]


"""Maggs-electrolyte methods"""


def get_electric_field_summary_stats(output_directory, temperature_index):
    """
    Returns the summary stats related to the zero mode of the electric field.

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
        The summary stats of the emergent electric field.  A one-dimensional numpy array of floats of length 4.
    """
    return np.loadtxt(f"{output_directory}/temp_{temperature_index:02d}/electric_field_summary_stats.csv",
                      dtype=float, delimiter=",")


def get_inverse_permittivity(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the inverse (electric) permittivity modulus
    [epsilon(x; temperature, no_of_sites)] ** (-1) =
    1.0 - 0.5 * no_of_sites * [\overline{E} - E[\overline{E}]] ** 2 / temperature, where E[.] is the expected value of
    the argument and \overline{E} is the zero mode of the electric field, as defined in Phys. Rev. B 91, 155412 (2015).

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_electric_field_summary_stats(output_directory, temperature_index)
    return 1.0 - summary_stats[0], summary_stats[1]


def get_topological_susceptibility(output_directory, temperature_index):
    r"""
    Returns the Monte Carlo mean and error of the topological susceptibility chi_w(x; temperature, no_of_sites) =
    0.5 * no_of_sites * [\overline{E}_w - E[\overline{E}_w]] ** 2 / temperature, where E[.] is the expected value of the
    argument and \overline{E}_w = 2 \pi w / L is the winding field, as defined in Phys. Rev. B 91, 155412 (2015).  In
    Phys. Rev. B 91, 155412 (2015), topological susceptibility was called the winding-field susceptibility and was
    defined without the factor of 1/2 to account for each Cartesian dimension of the system; in the return line below,
    this corresponds to the first np.mean() -> np.sum().

    Parameters
    ----------
    output_directory : str
        The location of the directory containing the sample(s) and Metropolis acceptance rate(s) (plurals refer to the
        option of multiple repeated simulations).
    temperature_index : int
        The index of the current sampling temperature within the configuration file.

    Returns
    -------
    Tuple[float, float]
        The tuple containing the mean and error.  A one-dimensional tuple of floats of length 2.  The 0th / 1st element
        is the Monte Carlo mean / error.
    """
    summary_stats = get_electric_field_summary_stats(output_directory, temperature_index)
    return summary_stats[2], summary_stats[3]
