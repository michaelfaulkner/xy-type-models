from version import version
import fileinput
import multiprocessing as mp
import os
import sys
import output.config_data_getter as config_data_getter


def main(executable, config_file):
    no_of_parallel_jobs = config_data_getter.get_basic_data(config_file)[8]
    if no_of_parallel_jobs < 1:
        print("ConfigurationError: For the value of no_of_parallel_jobs, give an integer not less than one.")
        exit()
    elif no_of_parallel_jobs == 1:
        print("Running a single Markov process.")
        run_single_simulation(executable, config_file)
    else:
        no_of_cpus = mp.cpu_count()
        if no_of_parallel_jobs < no_of_cpus:
            print("Running", no_of_parallel_jobs, "Markov processes in parallel on", no_of_parallel_jobs, "CPUs, where",
                  no_of_cpus, "CPUs are available.")
            pool = mp.Pool(no_of_parallel_jobs)
        else:
            print("Running", no_of_parallel_jobs, "Markov processes in parallel on", no_of_cpus, "CPUs, where",
                  no_of_cpus, "CPUs are available.")
            pool = mp.Pool(no_of_cpus)
        # create directory in which to store temporary copies of parent config file
        os.system("mkdir -p " + config_file.replace(".txt", ""))
        config_file_copies = [config_file.replace(".txt", "/job_" + str(job_number + 1) + ".txt") for job_number in
                              range(no_of_parallel_jobs)]
        for job_number, config_file_copy in enumerate(config_file_copies):
            # create temporary copies of parent config file
            os.system("cp " + config_file + " " + config_file_copy)
            for line in fileinput.input(config_file_copy, inplace=True):
                if 'output_directory' in line:
                    print(line.replace("' ", "/job_" + str(job_number + 1) + "'"), end="")
                else:
                    print(line, end="")
        pool.starmap(run_single_simulation, [(executable, config_file_copy) for config_file_copy in config_file_copies])
        pool.close()
        # delete temporary copies of parent config file
        os.system("rm -r " + config_file.replace(".txt", ""))


def print_start_message():
    """Print the start message."""
    print(f"xy-type-models (version {version}) - a Fortran/Python application for event-chain/Metropolis Monte Carlo "
          "simulation of the XY model, harmonic XY (HXY) model and Maggs lattice-field electrolyte model on a square, "
          "two-dimensional lattice (event-chain Monte Carlo only available for the XY and HXY models).'")


def run_single_simulation(executable, config_file):
    return os.system("./" + executable + " " + config_file)


if __name__ == "__main__":
    print_start_message()
    main(sys.argv[1], sys.argv[2])
