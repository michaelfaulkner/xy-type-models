import fileinput
import importlib
import multiprocessing as mp
import os
import sys


# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_modules = os.path.abspath(this_directory + '/output/')
sys.path.insert(0, directory_containing_modules)
config_data_getter = importlib.import_module('config_data_getter')


def main(executable, config_file_name):
    no_of_jobs = config_data_getter.get_basic_data(config_file_name)[8]
    if no_of_jobs < 1:
        print('ConfigurationError: For the value of no_of_jobs, give an integer not less than one.')
        exit()
    elif no_of_jobs == 1:
        run_single_simulation(executable, config_file_name)
    else:
        no_of_cpus = mp.cpu_count()
        if no_of_jobs < no_of_cpus:
            print("Running", no_of_jobs, "Markov processes in parallel on", no_of_jobs, "CPUs, where",
                  no_of_cpus, "CPUs are available.")
            pool = mp.Pool(no_of_jobs)
        else:
            print("Running", no_of_jobs, "Markov processes in parallel on", no_of_cpus, "CPUs, where",
                  no_of_cpus, "CPUs are available.")
            pool = mp.Pool(no_of_cpus)
        for job_number in range(no_of_jobs):
            os.system('cp ' + str(config_file_name) + ' ' + str(config_file_name) + '_run_' + str(job_number))
            for line in fileinput.input(str(config_file_name) + '_run_' + str(job_number)):
                if 'output_directory' in line:
                    print(line.replace("output_directory", "").replace(" ", "") + '_run_' + str(job_number))
                else:
                    print(line)
        config_file_copies = [str(config_file_name) + '_run_' + str(job_number) for job_number in no_of_jobs]
        pool.starmap(run_single_simulation, [(executable, config_file_copy) for config_file_copy in config_file_copies])
        pool.close()
        [os.system('rm -r ' + str(config_file_copy)) for config_file_copy in config_file_copies]


def run_single_simulation(executable, config_file_name):
    return os.system('./' + str(executable) + ' ' + str(config_file_name))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
