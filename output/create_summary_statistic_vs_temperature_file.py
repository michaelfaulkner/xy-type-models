import importlib
import multiprocessing as mp
import numpy as np
import os
import sys


# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file, summary_statistic_string):
    if (summary_statistic_string != "acceptance_rates" and summary_statistic_string != "no_of_events" and
            summary_statistic_string != "helicity_modulus" and summary_statistic_string != "magnetisation_norm" and
            summary_statistic_string != "magnetisation_phase" and summary_statistic_string != "specific_heat" and
            summary_statistic_string != "inverse_permittivity" and
            summary_statistic_string != "topological_sector_fluctuations" and
            summary_statistic_string != "inverse_vacuum_permittivity"):
        print("Give one of acceptance_rates, no_of_events, helicity_modulus, magnetisation_norm, magnetisation_phase, "
              "specific_heat, inverse_permittivity, topological_sector_fluctuations or inverse_vacuum_permittivity as "
              "the second positional argument.")
        exit()

    basic_config_data = config_data_getter.get_basic_data(config_file)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = (
        basic_config_data[0], basic_config_data[1], basic_config_data[2], basic_config_data[3], basic_config_data[5],
        basic_config_data[6], basic_config_data[7], basic_config_data[8])

    if (algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc") and summary_statistic_string == "acceptance_rates":
        print("ConfigurationError: These samples were generated by an event-chain Monte Carlo algorithm: rather than "
              "acceptance_rates, give no_of_events as the second positional argument.")
        exit()
    if ((algorithm_name == "xy-metropolis" or algorithm_name == "hxy-metropolis" or
         algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            summary_statistic_string == "no_of_events"):
        print("ConfigurationError: These samples were generated by a Metropolis algorithm: rather than "
              "no_of_events, give acceptance_rates as the second positional argument.")
        exit()
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (summary_statistic_string == "magnetisation_norm" or summary_statistic_string == "magnetisation_phase" or
             summary_statistic_string == "helicity_modulus" or
             summary_statistic_string == "inverse_vacuum_permittivity")):
        print("ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase, helicity_modulus or inverse_vacuum_permittivity as the second positional argument.")
        exit()
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-metropolis" or
         algorithm_name == "hxy-metropolis") and (summary_statistic_string == "inverse_permittivity" or
                                                  summary_statistic_string == "topological_sector_fluctuations")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity or "
              "topological_sector_fluctuations as the second positional argument.")
        exit()

    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / no_of_temperature_increments
    no_of_sites = integer_lattice_length ** 2

    temperature = initial_temperature
    output_file = open(output_directory + "/" + summary_statistic_string + "_vs_temperature.dat", "w")
    if summary_statistic_string == "acceptance_rates":
        temperature_directory = "/temp_eq_" + "{0:.2f}".format(temperature)
        acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        if len(acceptance_rates) == 2:
            output_file.write("temperature".ljust(15) + "final width of proposal interval".ljust(40) +
                              "rotational acceptance rate" + "\n")
        elif len(acceptance_rates) == 3:
            if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
                output_file.write("temperature".ljust(15) + "final width of proposal interval".ljust(40) +
                                  "acceptance rate (rotational moves)".ljust(40) +
                                  "acceptance rate (charge hops)" + "\n")
            else:
                output_file.write("temperature".ljust(15) + "final width of proposal interval".ljust(40) +
                                  "acceptance rate (rotational moves)".ljust(40) +
                                  "acceptance rate (external global moves)" + "\n")
        else:
            output_file.write("temperature".ljust(15) + "final width of proposal interval".ljust(40) +
                              "acceptance rate (field rotations)".ljust(40) +
                              "acceptance rate (charge hops)".ljust(40) + "acceptance rate (global moves)" + "\n")
    elif summary_statistic_string == "no_of_events":
        temperature_directory = "/temp_eq_" + "{0:.2f}".format(temperature)
        no_of_events = sample_getter.get_no_of_events(output_directory, temperature_directory)
        if len(no_of_events) == 1:
            output_file.write("temperature".ljust(15) + "number of events (field rotations)" + "\n")
        else:
            output_file.write("temperature".ljust(15) + "number of events (field rotations)".ljust(40) +
                              "acceptance rate (external global moves)" + "\n")
    else:
        output_file.write("temperature".ljust(15) + summary_statistic_string.ljust(35) + summary_statistic_string +
                          " error" + "\n")

    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = "/temp_eq_" + "{0:.2f}".format(temperature)
        if summary_statistic_string == "acceptance_rates" or summary_statistic_string == "no_of_events":
            get_sample_method = getattr(sample_getter, "get_" + summary_statistic_string)
            acceptance_rates_or_no_of_events = get_sample_method(output_directory, temperature_directory)
            if len(acceptance_rates_or_no_of_events) == 2:
                output_file.write("{0:.2f}".format(temperature).ljust(15) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[0]).ljust(40) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[1]) + "\n")
            elif len(acceptance_rates_or_no_of_events) == 3:
                output_file.write("{0:.2f}".format(temperature).ljust(15) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[0]).ljust(40) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[1]).ljust(40) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[2]) + "\n")
            else:
                output_file.write("{0:.2f}".format(temperature).ljust(15) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[0]).ljust(40) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[1]).ljust(40) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[2]).ljust(40) +
                                  "{0:.14e}".format(acceptance_rates_or_no_of_events[3]) + "\n")
        else:
            get_sample_method = getattr(sample_getter, "get_" + summary_statistic_string)
            if no_of_jobs == 1:
                sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)[
                         no_of_equilibration_sweeps:]
                sample_mean, sample_error = markov_chain_diagnostics.get_sample_mean_and_error(sample)
            else:
                number_of_cpus = mp.cpu_count()
                pool = mp.Pool(number_of_cpus)
                sample_means_and_errors = np.transpose(
                    np.array(pool.starmap(markov_chain_diagnostics.get_sample_mean_and_error,
                                          [[get_sample_method(output_directory + "/job_" + str(job_number + 1),
                                                              temperature_directory, beta, no_of_sites)[
                                            no_of_equilibration_sweeps:]] for job_number in range(no_of_jobs)])))
                pool.close()
                sample_mean = np.mean(sample_means_and_errors[0])
                sample_error = np.linalg.norm(sample_means_and_errors[1])
            output_file.write("{0:.2f}".format(temperature).ljust(15) + "{0:.14e}".format(sample_mean).ljust(35) +
                              "{0:.14e}".format(sample_error) + "\n")
        temperature += magnitude_of_temperature_increments
    output_file.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
