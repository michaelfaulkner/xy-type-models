import importlib
import matplotlib
import multiprocessing as mp
import os
import sys
import warnings

# import the run script as a module; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def check_for_observable_error(algorithm_name, observable_string):
    """Check to ensure that a valid observable has been given as the second positional argument."""
    if (observable_string != "potential" and observable_string != "specific_heat"
            and observable_string != "magnetisation_norm" and observable_string != "magnetisation_phase"
            and observable_string != "rotated_magnetisation_phase" and observable_string != "magnetisation_squared"
            and observable_string != "cartesian_magnetisation"
            and observable_string != "absolute_cartesian_magnetisation"
            and observable_string != "magnetic_susceptibility"
            and observable_string != "cartesian_relative_magnetisation"
            and observable_string != "relative_magnetisation_norm"
            and observable_string != "inverse_vacuum_permittivity" and observable_string != "macro_josephson_current"
            and observable_string != "helicity_modulus" and observable_string != "hxy_topological_sector"
            and observable_string != "hxy_internal_twist_susceptibility"
            and observable_string != "xy_twist_relaxation_susceptibility"
            and observable_string != "electric_field_zero_mode" and observable_string != "inverse_permittivity"
            and observable_string != "topological_sector" and observable_string != "topological_susceptibility"):
        print("ConfigurationError: Give one of potential, specific_heat, magnetisation_norm, magnetisation_phase, "
              "rotated_magnetisation_phase, magnetisation_squared, cartesian_magnetisation, "
              "absolute_cartesian_magnetisation, magnetic_susceptibility, cartesian_relative_magnetisation, "
              "relative_magnetisation_norm, inverse_vacuum_permittivity, macro_josephson_current, helicity_modulus, "
              "hxy_topological_sector, hxy_internal_twist_susceptibility, xy_twist_relaxation_susceptibility, "
              "electric_field_zero_mode, inverse_permittivity, topological_sector or topological_susceptibility as the "
              "second positional argument.")
        raise SystemExit
    check_for_observable_vs_model_error(algorithm_name, observable_string)


def check_for_observable_vs_model_error(algorithm_name, observable_string):
    """Raise an error if an XY or HXY observable has been given as the second positional argument for a
        Maggs-electrolyte model."""
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (observable_string == "magnetisation_norm" or observable_string == "magnetisation_phase"
             or observable_string == "rotated_magnetisation_phase" or observable_string == "magnetisation_squared"
             or observable_string == "cartesian_magnetisation"
             or observable_string == "absolute_cartesian_magnetisation"
             or observable_string == "magnetic_susceptibility"
             or observable_string == "cartesian_relative_magnetisation"
             or observable_string == "relative_magnetisation_norm" or observable_string == "inverse_vacuum_permittivity"
             or observable_string == "macro_josephson_current" or observable_string == "helicity_modulus"
             or observable_string == "hxy_topological_sector"
             or observable_string == "hxy_internal_twist_susceptibility"
             or observable_string == "xy_twist_relaxation_susceptibility")):
        print("ConfigurationError: This is a Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase, rotated_magnetisation_phase, magnetisation_squared, cartesian_magnetisation, "
              "absolute_cartesian_magnetisation, magnetic_susceptibility, cartesian_relative_magnetisation, "
              "relative_magnetisation_norm, inverse_vacuum_permittivity, macro_josephson_current, helicity_modulus, "
              "hxy_topological_sector, hxy_internal_twist_susceptibility or xy_twist_relaxation_susceptibility as the "
              "second positional argument.")
        raise SystemExit
    """Raise an error if a Maggs-electrolyte observable has been given as the second positional argument for an XY or 
        HXY model."""
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-uniform-noise-metropolis"
         or algorithm_name == "hxy-uniform-noise-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis"
         or algorithm_name == "hxy-gaussian-noise-metropolis") and (
            observable_string == "electric_field_zero_mode" or observable_string == "inverse_permittivity"
            or observable_string == "topological_sector" or observable_string == "topological_susceptibility")):
        print("ConfigurationError: This is an XY or HXY model: do not give either electric_field_zero_mode, "
              "inverse_permittivity, topological_sector or topological_susceptibility as the second positional "
              "argument.")
        raise SystemExit
    """Raise an error if hxy_topological_sector or hxy_internal_twist_susceptibility has been given as the second 
        positional argument for a non-HXY model."""
    if ((observable_string == "hxy_topological_sector" or observable_string == "hxy_internal_twist_susceptibility")
            and not (algorithm_name == "hxy-ecmc" or algorithm_name == "hxy-uniform-noise-metropolis"
                     or algorithm_name == "hxy-gaussian-noise-metropolis")):
        print("ConfigurationError: This is not the HXY model: do not give hxy_topological_sector or "
              "hxy_internal_twist_susceptibility as the second positional argument.")
        raise SystemExit
    """Raise an error if xy_twist_relaxation_field or xy_twist_relaxation_susceptibility has been given as the second 
            positional argument for a non-XY model."""
    if ((observable_string == "xy_twist_relaxation_field" or observable_string == "xy_twist_relaxation_susceptibility")
            and not (algorithm_name == "xy-ecmc" or algorithm_name == "xy-uniform-noise-metropolis"
                     or algorithm_name == "xy-gaussian-noise-metropolis")):
        print("ConfigurationError: This is not the XY model: do not give xy_twist_relaxation_field or "
              "xy_twist_relaxation_susceptibility as the second positional argument.")
        raise SystemExit


def setup_polyspectra_script(config_file, observable_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, _, temperatures, _,
     external_global_moves_string, no_of_runs, initial_run_index, max_no_of_cpus) = run_script.get_config_data(
        config_file)
    check_for_observable_error(algorithm_name, observable_string)
    check_initial_run_index(initial_run_index)
    return (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, temperatures,
            external_global_moves_string, no_of_runs, setup_pool(no_of_runs, max_no_of_cpus))


def setup_pool(no_of_runs, max_no_of_cpus):
    if no_of_runs < 1:
        raise Exception("ConfigurationError: For the value of no_of_runs, give an integer not less than one.")
    elif no_of_runs == 1:
        print("Running a single sample-analysis process.")
        return None
    else:
        no_of_available_cpus = mp.cpu_count()
        if no_of_available_cpus > max_no_of_cpus:
            no_of_cpus = max_no_of_cpus
        else:
            no_of_cpus = no_of_available_cpus
        if no_of_runs < no_of_cpus:
            print(f"Running {no_of_runs} sample-analysis processes in parallel on {no_of_runs} CPUs, "
                  f"where {no_of_available_cpus} CPUs are available.")
            pool = mp.Pool(no_of_runs)
        else:
            print(f"Running {no_of_runs} sample-analysis processes in parallel on {no_of_cpus} CPUs, where "
                  f"{no_of_available_cpus} CPUs are available.")
            pool = mp.Pool(no_of_cpus)
        return pool


def get_sample_is_one_dimensional(observable_string):
    if (observable_string == "potential" or observable_string == "specific_heat"
            or observable_string == "magnetisation_norm" or observable_string == "magnetisation_phase"
            or observable_string == "rotated_magnetisation_phase" or observable_string == "magnetisation_squared"
            or observable_string == "magnetic_susceptibility" or observable_string == "relative_magnetisation_norm"
            or observable_string == "inverse_vacuum_permittivity" or observable_string == "helicity_modulus"
            or observable_string == "hxy_internal_twist_susceptibility"
            or observable_string == "xy_twist_relaxation_susceptibility" or observable_string == "inverse_permittivity"
            or observable_string == "topological_susceptibility"):
        return True
    else:
        return False


def reverse_enumerate(iterable_object):
    return [(len(iterable_object) - index - 1, element) for index, element in
            enumerate(list(reversed(iterable_object)))]


def check_initial_run_index(initial_run_index):
    if initial_run_index != 0:
        warnings.warn("Warning: The value of initial_run_index is not equal to zero.  Configuration files of this type "
                      "are designed for running Fortran code only.")
