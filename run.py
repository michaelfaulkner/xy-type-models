from version import version
import fileinput
import multiprocessing as mp
import os
import sys
import output.config_data_getter as config_data_getter


def main(config_file):
    basic_config_data = config_data_getter.get_basic_data(config_file)
    algorithm_name, no_of_parallel_jobs = basic_config_data[0], basic_config_data[8]

    if algorithm_name == "hxy-ecmc":
        executable = "executables/hxy_ecmc_algorithm.exe"
    elif algorithm_name == "hxy-metropolis":
        executable = "executables/hxy_metropolis_algorithm.exe"
    elif algorithm_name == "xy-ecmc":
        executable = "executables/xy_ecmc_algorithm.exe"
    elif algorithm_name == "xy-metropolis":
        executable = "executables/xy_metropolis_algorithm.exe"
    elif algorithm_name == "elementary-electrolyte":
        executable = "executables/elementary_electrolyte_algorithm.exe"
    elif algorithm_name == "multivalued-electrolyte":
        executable = "executables/multivalued_electrolyte_algorithm.exe"
    else:
        executable = None
        print("ConfigurationError: For the value of algorithm_name, give one of hxy-ecmc, hxy-metropolis, xy-ecmc, "
              "xy-metropolis, elementary-electrolyte or multivalued-electrolyte.")
        exit()

    if not os.path.isfile(executable):
        print(f"Executable for the {algorithm_name} algorithm does not exist. Run 'make' or 'make {algorithm_name}' "
              f"then try again.")
        exit()

    if no_of_parallel_jobs < 1:
        print("ConfigurationError: For the value of no_of_parallel_jobs, give an integer not less than one.")
        exit()
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
        # create directory in which to store temporary copies of parent config file
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
    print(f"xy-type-models (version {version}) - https://github.com/michaelfaulkner/xy-type-models/ - a Fortran/Python "
          f"application for event-chain/Metropolis Monte Carlo simulation of the XY model, harmonic XY (HXY) model and "
          f"Maggs lattice-field electrolyte model on a square, two-dimensional lattice (event-chain Monte Carlo only "
          f"available for the XY and HXY models).'")


def run_single_simulation(executable, config_file):
    return os.system(f"./{executable} {config_file}")


if __name__ == "__main__":
    print_start_message()
    main(sys.argv[1])
