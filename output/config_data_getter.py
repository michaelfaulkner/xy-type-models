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


def check_for_config_errors(algorithm_name, observable_string):
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (observable_string == "magnetisation_norm" or observable_string == "magnetisation_norm" or
             observable_string == "helicity_modulus" or observable_string == "inverse_vacuum_permittivity" or
             observable_string == "toroidal_vortex_polarisation")):
        print("ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase, helicity_modulus, inverse_vacuum_permittivity or toroidal_vortex_polarisation as "
              "the second positional argument.")
        raise SystemExit
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-metropolis" or
         algorithm_name == "hxy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis" or
         algorithm_name == "hxy-gaussian-noise-metropolis") and (
            observable_string == "inverse_permittivity" or
            observable_string == "topological_sector_fluctuations" or
            observable_string == "toroidal_polarisation")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity, "
              "topological_sector_fluctuations or toroidal_polarisation as the second positional argument.")
        raise SystemExit


def get_current_temperature_and_magnitude_of_increments(final_temperature, initial_temperature,
                                                        no_of_temperature_increments):
    temperature = final_temperature
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / no_of_temperature_increments
    return temperature, magnitude_of_temperature_increments


def get_polyspectra_config_data(config_file, observable_string):
    basic_config_data = get_basic_data(config_file)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = (basic_config_data[0], basic_config_data[1],
                                                                     basic_config_data[2], basic_config_data[3],
                                                                     basic_config_data[5], basic_config_data[6],
                                                                     basic_config_data[7], basic_config_data[8])
    check_for_config_errors(algorithm_name, observable_string)
    no_of_sites = integer_lattice_length ** 2
    return (final_temperature, initial_temperature, no_of_sites, no_of_equilibration_sweeps, no_of_jobs,
            no_of_temperature_increments, output_directory)
