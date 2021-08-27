from version import version
import fileinput
import multiprocessing as mp
import os
import sys
import output.setup_scripts as setup_scripts


def main(config_file):
    config_data = setup_scripts.get_config_data(config_file)
    algorithm_name, no_of_parallel_jobs = config_data[0], config_data[7]
    executable = get_executable(algorithm_name)
    if no_of_parallel_jobs < 1:
        raise Exception("ConfigurationError: For the value of no_of_parallel_jobs, give an integer not less than one.")
    elif no_of_parallel_jobs == 1:
        print("Running a single Markov process.")
        run_single_simulation(executable, config_file)
    else:
        no_of_cpus = mp.cpu_count()
        if no_of_parallel_jobs < no_of_cpus:
            print(f"Running {no_of_parallel_jobs} Markov processes in parallel on {no_of_parallel_jobs} CPUs, where",
                  f"{no_of_cpus} CPUs are available.")
            pool = mp.Pool(no_of_parallel_jobs)
        else:
            print(f"Running {no_of_parallel_jobs} Markov processes in parallel on {no_of_cpus} CPUs, where "
                  f"{no_of_cpus} CPUs are available.")
            pool = mp.Pool(no_of_cpus)
        # create directory in which to store temporary copies of the parent config file
        os.system(f"mkdir -p {config_file.replace('.txt', '')}")
        config_file_copies = [config_file.replace(".txt", f"/job_{job_number + 1}.txt") for job_number in
                              range(no_of_parallel_jobs)]
        for job_number, config_file_copy in enumerate(config_file_copies):
            # create temporary copies of parent config file
            os.system(f"cp {config_file} {config_file_copy}")
            for line in fileinput.input(config_file_copy, inplace=True):
                if "output_directory" in line:
                    print(line.replace("' ", f"/job_{job_number + 1}'"), end="")
                else:
                    print(line, end="")
        pool.starmap(run_single_simulation, [(executable, config_file_copy) for config_file_copy in config_file_copies])
        pool.close()
        # delete temporary copies of parent config file
        os.system(f"rm -r {config_file.replace('.txt', '')}")


def print_start_message():
    """Print the start message."""
    print(f"xy-type-models (version {version}) - https://github.com/michaelfaulkner/xy-type-models/ - a hybrid "
          f"Fortran-Python application for event-chain/Metropolis Monte Carlo simulation of two-dimensional XY-type "
          f"models in statistical physics.")


def get_executable(algorithm_name):
    if algorithm_name == "hxy-ecmc":
        executable = "executables/hxy_ecmc_algorithm.exe"
    elif algorithm_name == "hxy-metropolis":
        executable = "executables/hxy_metropolis_algorithm.exe"
    elif algorithm_name == "hxy-gaussian-noise-metropolis":
        executable = "executables/hxy_gaussian_noise_metropolis_algorithm.exe"
    elif algorithm_name == "xy-ecmc":
        executable = "executables/xy_ecmc_algorithm.exe"
    elif algorithm_name == "xy-metropolis":
        executable = "executables/xy_metropolis_algorithm.exe"
    elif algorithm_name == "xy-gaussian-noise-metropolis":
        executable = "executables/xy_gaussian_noise_metropolis_algorithm.exe"
    elif algorithm_name == "elementary-electrolyte":
        executable = "executables/elementary_electrolyte_algorithm.exe"
    elif algorithm_name == "multivalued-electrolyte":
        executable = "executables/multivalued_electrolyte_algorithm.exe"
    else:
        raise Exception("ConfigurationError: For the value of algorithm_name, give one of hxy-ecmc, hxy-metropolis, "
                        "xy-ecmc, xy-metropolis, elementary-electrolyte or multivalued-electrolyte.")
    if not os.path.isfile(executable):
        raise Exception(f"Executable for the {algorithm_name} algorithm does not exist.  Run 'make' or 'make "
                        f"{algorithm_name}' then try again.")
    return executable


def run_single_simulation(executable, config_file):
    return os.system(f"./{executable} {config_file}")


if __name__ == "__main__":
    print_start_message()
    main(sys.argv[1])
