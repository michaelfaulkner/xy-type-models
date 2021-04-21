import csv


def get_basic_configuration_data(config_file_name):
    with open(config_file_name, 'r') as config_file:
        config_data = csv.reader(config_file, delimiter='\t')
        for index, row in enumerate(config_data):
            if index == 0:
                algorithm_name = row[0].replace("'", "").replace("algorithm_name", "").replace(" ", "")
            if index == 1:
                simulation_directory = row[0].replace("'", "").replace("output directory", "").replace(" ", "")
            if index == 2:
                lattice_length = int(row[0].replace("'", "").replace("lattice length", "").replace(" ", ""))
            if index == 3:
                no_of_equilibrium_iterations = int(row[0].replace("number of equilibration iterations", "").replace(" ",
                                                                                                                    ""))
            if index == 4:
                no_of_observations = int(row[0].replace("number of observations (sample size)", "").replace(" ", ""))
            if index == 5:
                initial_temperature = float(row[0].replace("initial temperature", "").replace(" ", ""))
            if index == 6:
                final_temperature = float(row[0].replace("final temperature", "").replace(" ", ""))
            if index == 7:
                no_of_temperature_increments = int(row[0].replace("number of temperature increments", "").replace(" ",
                                                                                                                  ""))
                break
    return (algorithm_name, simulation_directory, lattice_length, no_of_equilibrium_iterations, no_of_observations,
            initial_temperature, final_temperature, no_of_temperature_increments)
