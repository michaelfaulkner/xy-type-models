import csv
import matplotlib
import multiprocessing as mp


def get_config_data(config_file_location):
    with open(config_file_location, 'r') as config_file:
        for row in csv.reader(config_file, delimiter='\t'):
            if 'algorithm_name' in row[0]:
                algorithm_name = row[0].replace("'", "").replace("algorithm_name", "").replace(" ", "")
            if 'output_directory' in row[0]:
                output_directory = row[0].replace("'", "").replace("output_directory", "").replace(" ", "")
            if 'integer_lattice_length' in row[0]:
                integer_lattice_length = int(row[0].replace("'", "").replace("integer_lattice_length", "").replace(" ",
                                                                                                                   ""))
                no_of_sites = integer_lattice_length ** 2
            if 'no_of_equilibration_sweeps' in row[0]:
                no_of_equilibration_sweeps = int(row[0].replace("no_of_equilibration_sweeps", "").replace(" ", ""))
            if 'initial_temperature' in row[0]:
                initial_temperature = float(row[0].replace("d0", "").replace("initial_temperature", "").replace(" ",
                                                                                                                ""))
            if 'final_temperature' in row[0]:
                final_temperature = float(row[0].replace("d0", "").replace("final_temperature", "").replace(" ", ""))
            if 'no_of_temperature_increments' in row[0]:
                no_of_temperature_increments = int(row[0].replace("no_of_temperature_increments", "").replace(" ", ""))
            if 'no_of_parallel_jobs' in row[0]:
                no_of_parallel_jobs = int(row[0].replace("no_of_parallel_jobs", "").replace(" ", ""))
    return (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
            final_temperature, no_of_temperature_increments, no_of_parallel_jobs)


def check_for_observable_error(algorithm_name, observable_string):
    if (observable_string != "magnetisation_norm" and observable_string != "magnetisation_phase" and
            observable_string != "cartesian_magnetisation" and observable_string != "helicity_modulus" and
            observable_string != "inverse_vacuum_permittivity" and observable_string != "toroidal_vortex_polarisation"
            and observable_string != "specific_heat" and observable_string != "inverse_permittivity" and
            observable_string != "topological_sector_fluctuations" and
            observable_string != "toroidal_vortex_polarisation"):
        print("ConfigurationError: Give one of magnetisation_norm, magnetisation_phase, cartesian_magnetisation, "
              "helicity_modulus, inverse_vacuum_permittivity, toroidal_vortex_polarisation, specific_heat, "
              "inverse_permittivity, topological_sector_fluctuations or toroidal_vortex_polarisation as the second "
              "positional argument.")
        raise SystemExit
    check_for_observable_vs_algorithm_error(algorithm_name, observable_string)


def check_for_observable_vs_algorithm_error(algorithm_name, observable_string):
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (observable_string == "magnetisation_norm" or observable_string == "magnetisation_phase" or
             observable_string == "cartesian_magnetisation" or observable_string == "helicity_modulus" or
             observable_string == "inverse_vacuum_permittivity" or
             observable_string == "toroidal_vortex_polarisation")):
        print("ConfigurationError: This is a Maggs-electrolyte model: do not give either magnetisation_norm, "
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


def get_temperature_and_magnitude_of_increments(initial_temperature, final_temperature, no_of_temperature_increments):
    temperature = final_temperature
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
    return temperature, magnitude_of_temperature_increments


def set_up_polyspectra_script(config_file, observable_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = get_config_data(config_file)
    check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)
    if no_of_jobs > 1:
        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)
    else:
        pool = None
    return (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
            no_of_jobs, temperature, magnitude_of_temperature_increments, pool)