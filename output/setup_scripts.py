import importlib
import matplotlib
import multiprocessing as mp
import os
import sys

# import the run script as a module; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def check_for_observable_error(algorithm_name, observable_string):
    if (observable_string != "magnetisation_norm" and observable_string != "magnetisation_phase"
            and observable_string != "cartesian_magnetisation" and observable_string != "magnetic_susceptibility"
            and observable_string != "helicity_modulus" and observable_string != "inverse_vacuum_permittivity"
            and observable_string != "toroidal_vortex_polarisation" and observable_string != "specific_heat"
            and observable_string != "inverse_permittivity" and observable_string != "topological_sector_fluctuations"
            and observable_string != "toroidal_vortex_polarisation"):
        print("ConfigurationError: Give one of magnetisation_norm, magnetisation_phase, cartesian_magnetisation, "
              "magnetic_susceptibility, helicity_modulus, inverse_vacuum_permittivity, toroidal_vortex_polarisation, "
              "specific_heat, inverse_permittivity, topological_sector_fluctuations or toroidal_vortex_polarisation as "
              "the second positional argument.")
        raise SystemExit
    check_for_observable_vs_algorithm_error(algorithm_name, observable_string)


def check_for_observable_vs_algorithm_error(algorithm_name, observable_string):
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (observable_string == "magnetisation_norm" or observable_string == "magnetisation_phase"
             or observable_string == "cartesian_magnetisation" or observable_string == "magnetic_susceptibility"
             or observable_string == "helicity_modulus" or observable_string == "inverse_vacuum_permittivity"
             or observable_string == "toroidal_vortex_polarisation")):
        print("ConfigurationError: This is a Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase, magnetic_susceptibility, helicity_modulus, inverse_vacuum_permittivity or "
              "toroidal_vortex_polarisation as the second positional argument.")
        raise SystemExit
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-metropolis" or
         algorithm_name == "hxy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis" or
         algorithm_name == "hxy-gaussian-noise-metropolis") and (
            observable_string == "inverse_permittivity" or observable_string == "topological_sector_fluctuations" or
            observable_string == "toroidal_polarisation")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity, "
              "topological_sector_fluctuations or toroidal_polarisation as the second positional argument.")
        raise SystemExit


def get_temperature_and_magnitude_of_increments(initial_temperature, final_temperature, no_of_temperature_increments):
    temperature = final_temperature
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
    return temperature, magnitude_of_temperature_increments


def set_up_polyspectra_script(config_file, observable_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature, final_temperature,
     no_of_temperature_increments, no_of_jobs, max_no_of_cpus) = run_script.get_config_data(config_file)
    check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)
    if no_of_jobs < 1:
        raise Exception("ConfigurationError: For the value of no_of_jobs, give an integer not less than one.")
    elif no_of_jobs == 1:
        print("Running a single sample-analysis process.")
        pool = None
    else:
        no_of_available_cpus = mp.cpu_count()
        if no_of_available_cpus > max_no_of_cpus:
            no_of_cpus = max_no_of_cpus
        else:
            no_of_cpus = no_of_available_cpus
        if no_of_jobs < no_of_cpus:
            print(f"Running {no_of_jobs} sample-analysis processes in parallel on {no_of_jobs} CPUs, "
                  f"where {no_of_available_cpus} CPUs are available.")
            pool = mp.Pool(no_of_jobs)
        else:
            print(f"Running {no_of_jobs} sample-analysis processes in parallel on {no_of_cpus} CPUs, where "
                  f"{no_of_available_cpus} CPUs are available.")
            pool = mp.Pool(no_of_cpus)
    return (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
            no_of_jobs, temperature, magnitude_of_temperature_increments, pool)
