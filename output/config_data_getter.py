import csv


def get_basic_data(config_file_location):
    with open(config_file_location, 'r') as config_file:
        for row in csv.reader(config_file, delimiter='\t'):
            if 'algorithm_name' in row[0]:
                algorithm_name = row[0].replace("'", "").replace("algorithm_name", "").replace(" ", "")
            if 'output_directory' in row[0]:
                output_directory = row[0].replace("'", "").replace("output_directory", "").replace(" ", "")
            if 'integer_lattice_length' in row[0]:
                integer_lattice_length = int(row[0].replace("'", "").replace("integer_lattice_length", "").replace(" ",
                                                                                                                   ""))
            if 'no_of_equilibration_sweeps' in row[0]:
                no_of_equilibration_sweeps = int(row[0].replace("no_of_equilibration_sweeps", "").replace(" ", ""))
            if 'no_of_observations' in row[0]:
                no_of_observations = int(row[0].replace("no_of_observations", "").replace(" ", ""))
            if 'initial_temperature' in row[0]:
                initial_temperature = float(row[0].replace("d0", "").replace("initial_temperature", "").replace(" ",
                                                                                                                ""))
            if 'final_temperature' in row[0]:
                final_temperature = float(row[0].replace("d0", "").replace("final_temperature", "").replace(" ", ""))
            if 'no_of_temperature_increments' in row[0]:
                no_of_temperature_increments = int(row[0].replace("no_of_temperature_increments", "").replace(" ", ""))
            if 'no_of_parallel_jobs' in row[0]:
                no_of_parallel_jobs = int(row[0].replace("no_of_parallel_jobs", "").replace(" ", ""))
    return (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, no_of_observations,
            initial_temperature, final_temperature, no_of_temperature_increments, no_of_parallel_jobs)
