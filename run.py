"""Executable script which runs the xy-type-models application based on an input configuration file."""
from version import version
import csv
import fileinput
import multiprocessing as mp
import os
import random
import sys


def print_start_message():
    """Print the start message."""
    print(f"xy-type-models (version {version}) - https://github.com/michaelfaulkner/xy-type-models/ - a "
          f"Fortran-Python application for event-chain/Metropolis Monte Carlo simulation of two-dimensional XY-type "
          f"models in statistical physics.")


def main(config_file_location):
    """
    Use the location of the configuration file to run the xy-type-models application.

    Parameters
    ----------
    config_file_location : str
        The location of the configuration file.

    Returns
    -------
    None
    """
    config_data = get_config_data(config_file_location)
    algorithm_name, no_of_runs, initial_run_index, max_no_of_cpus = (config_data[0], config_data[9], config_data[10],
                                                                     config_data[11])
    executable_location = get_executable(algorithm_name)
    if no_of_runs < 1:
        raise Exception("ConfigurationError: For the value of no_of_runs, give an integer not less than one.")
    elif no_of_runs == 1:
        print("Running a single Markov process.")
        random_seed = random.randint(100000000, 999999999)
        run_single_simulation(executable_location, random_seed, config_file_location)
    else:
        no_of_available_cpus = mp.cpu_count()
        if no_of_available_cpus > max_no_of_cpus:
            no_of_cpus = max_no_of_cpus
        else:
            no_of_cpus = no_of_available_cpus
        if no_of_runs < no_of_cpus:
            print(f"Running {no_of_runs} Markov processes in parallel on {no_of_runs} CPUs, where",
                  f"{no_of_available_cpus} CPUs are available.")
            pool = mp.Pool(no_of_runs)
        else:
            print(f"Running {no_of_runs} Markov processes in parallel on {no_of_cpus} CPUs, where "
                  f"{no_of_available_cpus} CPUs are available.")
            pool = mp.Pool(no_of_cpus)
        # create directory in which to store temporary copies of the parent config file
        os.system(f"mkdir -p {config_file_location.replace('.txt', '')}")
        random_seeds = [random.randint(100000000, 999999999) for _ in range(no_of_runs)]
        config_file_copies = [config_file_location.replace(".txt", f"/run_{run_number + initial_run_index}.txt")
                              for run_number in range(no_of_runs)]
        for run_number, config_file_copy in enumerate(config_file_copies):
            # create temporary copies of parent config file
            os.system(f"cp {config_file_location} {config_file_copy}")
            for line in fileinput.input(config_file_copy, inplace=True):
                if "output_directory" in line:
                    print(line.replace("' ", f"/run_{run_number + initial_run_index}'"), end="")
                else:
                    print(line, end="")
        pool.starmap(run_single_simulation, [(executable_location, random_seeds[run_number], config_file_copy)
                                             for run_number, config_file_copy in enumerate(config_file_copies)])
        pool.close()
        # delete temporary copies of parent config file
        os.system(f"rm -r {config_file_location.replace('.txt', '')}")


def get_config_data(config_file_location):
    """
    Returns the data from within the configuration file.

    Parameters
    ----------
    config_file_location : str
        The location of the configuration file.

    Returns
    -------
    (str, str, int, int, float, float, int, int, int)
        The configuration-file data.  A one-dimensional tuple of length 9.
    """
    with open(config_file_location, 'r') as config_file:
        for row in csv.reader(config_file, delimiter='\t'):
            if 'algorithm_name' in row[0]:
                algorithm_name = row[0].replace("'", "").replace("algorithm_name", "").replace(" ", "")
            if 'output_directory' in row[0]:
                output_directory = row[0].replace("'", "").replace("output_directory", "").replace(" ", "")
            if 'integer_lattice_length' in row[0]:
                integer_lattice_length = int(row[0].replace("'", "").replace("integer_lattice_length", "").replace(" ",
                                                                                                                   ""))
                if algorithm_name == "3dxy-gaussian-noise-metropolis":
                    no_of_sites = integer_lattice_length ** 3
                    no_of_sites_string = (f"{integer_lattice_length}x{integer_lattice_length}x{integer_lattice_length}_"
                                          f"sites")
                else:
                    no_of_sites = integer_lattice_length ** 2
                    no_of_sites_string = f"{integer_lattice_length}x{integer_lattice_length}_sites"
            if 'no_of_equilibration_sweeps' in row[0]:
                no_of_equilibration_sweeps = int(row[0].replace("no_of_equilibration_sweeps", "").replace(" ", ""))
            if 'no_of_samples' in row[0]:
                no_of_samples = int(row[0].replace("no_of_samples", "").replace(" ", ""))
            if 'initial_temperature' in row[0]:
                initial_temperature = float(row[0].replace("d0", "").replace("initial_temperature", "").replace(" ",
                                                                                                                ""))
            if 'final_temperature' in row[0]:
                final_temperature = float(row[0].replace("d0", "").replace("final_temperature", "").replace(" ", ""))
            if 'no_of_temperature_increments' in row[0]:
                no_of_temperature_increments = int(row[0].replace("no_of_temperature_increments", "").replace(" ", ""))
                if type(no_of_temperature_increments) is not int or not (0 <= no_of_temperature_increments <= 99):
                    raise Exception("ConfigurationError: The value of no_of_temperature_increments must be an integer "
                                    "not less than zero and not greater than 99.")
            if 'use_external_global_moves' in row[0]:
                if '.true.' in row[0]:
                    use_external_global_moves = True
                    if 'electrolyte' in algorithm_name:
                        external_global_moves_string = "w_global_moves"
                    else:
                        external_global_moves_string = "w_twists"
                else:
                    use_external_global_moves = False
                    if 'electrolyte' in algorithm_name:
                        external_global_moves_string = "sans_global_moves"
                    else:
                        external_global_moves_string = "sans_twists"
            if 'no_of_runs' in row[0]:
                no_of_runs = int(row[0].replace("no_of_runs", "").replace(" ", ""))
            if 'initial_run_index' in row[0]:
                initial_run_index = int(row[0].replace("initial_run_index", "").replace(" ", ""))
            if 'max_no_of_cpus' in row[0]:
                max_no_of_cpus = int(row[0].replace("max_no_of_cpus", "").replace(" ", ""))
    return (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
            no_of_samples, get_temperatures(initial_temperature, final_temperature, no_of_temperature_increments),
            use_external_global_moves, external_global_moves_string, no_of_runs, initial_run_index, max_no_of_cpus)


def get_temperatures(initial_temperature, final_temperature, no_of_temperature_increments):
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
    return [initial_temperature + temperature_index * magnitude_of_temperature_increments
            for temperature_index in range(no_of_temperature_increments + 1)]


def get_executable(algorithm_name):
    """
    Returns the location of the (Fortran) sampling-algorithm executable.

    Parameters
    ----------
    algorithm_name : str
        The name of the required sampling algorithm.

    Returns
    -------
    str
        The location of the sampling-algorithm executable.
    """
    if algorithm_name == "hxy-ecmc":
        executable_location = "executables/hxy_ecmc_algorithm.exe"
    elif algorithm_name == "hxy-uniform-noise-metropolis":
        executable_location = "executables/hxy_uniform_noise_metropolis_algorithm.exe"
    elif algorithm_name == "hxy-gaussian-noise-metropolis":
        executable_location = "executables/hxy_gaussian_noise_metropolis_algorithm.exe"
    elif algorithm_name == "xy-ecmc":
        executable_location = "executables/xy_ecmc_algorithm.exe"
    elif algorithm_name == "xy-uniform-noise-metropolis":
        executable_location = "executables/xy_uniform_noise_metropolis_algorithm.exe"
    elif algorithm_name == "xy-gaussian-noise-metropolis":
        executable_location = "executables/xy_gaussian_noise_metropolis_algorithm.exe"
    elif algorithm_name == "3dxy-gaussian-noise-metropolis":
        executable_location = "executables/3dxy_gaussian_noise_metropolis_algorithm.exe"
    elif algorithm_name == "elementary-electrolyte":
        executable_location = "executables/elementary_electrolyte_algorithm.exe"
    elif algorithm_name == "multivalued-electrolyte":
        executable_location = "executables/multivalued_electrolyte_algorithm.exe"
    else:
        raise Exception("ConfigurationError: For the value of algorithm_name, give one of hxy-ecmc, "
                        "hxy-uniform-noise-metropolis, hxy-gaussian-noise-metropolis, xy-ecmc, "
                        "xy-uniform-noise-metropolis, xy-gaussian-noise-metropolis, 3dxy-gaussian-noise-metropolis, "
                        "elementary-electrolyte or multivalued-electrolyte.")
    if not os.path.isfile(executable_location):
        raise Exception(f"Executable for the {algorithm_name} algorithm does not exist.  Run 'make' or 'make "
                        f"{algorithm_name}' then try again.")
    return executable_location


def run_single_simulation(executable_location, random_seed, config_file_location):
    """
    Runs a single simulation using the (Fortran) sampling algorithm.

    Parameters
    ----------
    executable_location : str
        The location of the Fortran sampling-algorithm executable.
    random_seed : int
        The random seed for the random-number generator in the Fortran executable.  A nine-digit integer.
    config_file_location : str
        The location of the configuration file.

    Returns
    -------
    None
    """
    os.system(f"./{executable_location} {random_seed} {config_file_location}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("InterfaceError: One positional argument required - give the location of the configuration "
                        "file.")
    print_start_message()
    main(sys.argv[1])
