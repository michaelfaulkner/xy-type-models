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
    basic_config_data = config_data_getter.get_basic_data(config_file)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = (
        basic_config_data[0], basic_config_data[1], basic_config_data[2], basic_config_data[3], basic_config_data[5],
        basic_config_data[6], basic_config_data[7], basic_config_data[8])

    check_for_config_errors(algorithm_name, summary_statistic_string)
    temperature = initial_temperature
    no_of_sites = integer_lattice_length ** 2
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / no_of_temperature_increments

    output_file = open(output_directory + "/" + summary_statistic_string + "_vs_temperature.dat", "w")
    if summary_statistic_string == "acceptance_rates":
        temperature_directory = f"temp_eq_{temperature:.2f}"
        if no_of_jobs == 1:
            acceptance_rates = sample_getter.get_acceptance_rates(output_directory, temperature_directory)
        else:
            acceptance_rates = sample_getter.get_acceptance_rates(output_directory + "/job_1", temperature_directory)
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
        temperature_directory = f"temp_eq_{temperature:.2f}"
        if no_of_jobs == 1:
            no_of_events = sample_getter.get_no_of_events(output_directory, temperature_directory)
        else:
            no_of_events = sample_getter.get_no_of_events(output_directory + "/job_1", temperature_directory)
        if len(no_of_events) == 1:
            output_file.write("temperature".ljust(15) + "number of events (field rotations)" + "\n")
        else:
            output_file.write("temperature".ljust(15) + "number of events (field rotations)".ljust(40) +
                              "acceptance rate (external global moves)" + "\n")
    else:
        output_file.write("temperature".ljust(15) + summary_statistic_string.ljust(35) + summary_statistic_string +
                          " error" + "\n")

    if no_of_jobs > 1:
        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)

    for i in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"
        if summary_statistic_string == "acceptance_rates" or summary_statistic_string == "no_of_events":
            get_sample_method = getattr(sample_getter, "get_" + summary_statistic_string)
            if no_of_jobs == 1:
                acceptance_rates_or_no_of_events = get_sample_method(output_directory, temperature_directory)
            else:
                acceptance_rates_or_no_of_events = np.mean([
                    get_sample_method(output_directory + "/job_" + str(job_number + 1), temperature_directory)
                    for job_number in range(no_of_jobs)], axis=0)
            if len(acceptance_rates_or_no_of_events) == 1:
                output_file.write(f"{temperature:.2f}".ljust(15) +
                                  f"{acceptance_rates_or_no_of_events[0]:.14e}" + "\n")
            elif len(acceptance_rates_or_no_of_events) == 2:
                output_file.write(f"{temperature:.2f}".ljust(15) +
                                  f"{acceptance_rates_or_no_of_events[0]:.14e}".ljust(40) +
                                  f"{acceptance_rates_or_no_of_events[1]:.14e}" + "\n")
            elif len(acceptance_rates_or_no_of_events) == 3:
                output_file.write(f"{temperature:.2f}".ljust(15) +
                                  f"{acceptance_rates_or_no_of_events[0]:.14e}".ljust(40) +
                                  f"{acceptance_rates_or_no_of_events[1]:.14e}".ljust(40) +
                                  f"{acceptance_rates_or_no_of_events[2]:.14e}" + "\n")
            else:
                output_file.write(f"{temperature:.2f}".ljust(15) +
                                  f"{acceptance_rates_or_no_of_events[0]:.14e}".ljust(40) +
                                  f"{acceptance_rates_or_no_of_events[1]:.14e}".ljust(40) +
                                  f"{acceptance_rates_or_no_of_events[2]:.14e}".ljust(40) +
                                  f"{acceptance_rates_or_no_of_events[3]:.14e}" + "\n")
        else:
            get_sample_method = getattr(sample_getter, "get_" + summary_statistic_string)
            if no_of_jobs == 1:
                sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)[
                         no_of_equilibration_sweeps:]
                sample_mean, sample_error = markov_chain_diagnostics.get_sample_mean_and_error(sample)
            else:
                sample_means_and_errors = np.transpose(
                    np.array(pool.starmap(markov_chain_diagnostics.get_sample_mean_and_error,
                                          [[get_sample_method(output_directory + "/job_" + str(job_number + 1),
                                                              temperature_directory, beta, no_of_sites)[
                                            no_of_equilibration_sweeps:]] for job_number in range(no_of_jobs)])))
                sample_mean = np.mean(sample_means_and_errors[0])
                sample_error = np.linalg.norm(sample_means_and_errors[1])
            output_file.write(f"{temperature:.2f}".ljust(15) + f"{sample_mean:.14e}".ljust(35) + f"{sample_error:.14e}"
                              + "\n")
        temperature += magnitude_of_temperature_increments
    if no_of_jobs > 1:
        pool.close()
    output_file.close()


def check_for_config_errors(algorithm_name, summary_statistic_string):
    if (summary_statistic_string != "acceptance_rates" and summary_statistic_string != "no_of_events" and
            summary_statistic_string != "helicity_modulus" and summary_statistic_string != "magnetisation_norm" and
            summary_statistic_string != "magnetisation_phase" and summary_statistic_string != "specific_heat" and
            summary_statistic_string != "inverse_permittivity" and
            summary_statistic_string != "topological_sector_fluctuations" and
            summary_statistic_string != "inverse_vacuum_permittivity"):
        print("ConfigurationError: Give one of acceptance_rates, no_of_events, helicity_modulus, magnetisation_norm, "
              "magnetisation_phase, specific_heat, inverse_permittivity, topological_sector_fluctuations or "
              "inverse_vacuum_permittivity as the second positional argument.")
        exit()
    if (algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc") and summary_statistic_string == "acceptance_rates":
        print("ConfigurationError: These samples were generated by an event-chain Monte Carlo algorithm: rather than "
              "acceptance_rates, give no_of_events as the second positional argument.")
        exit()
    if ((algorithm_name == "xy-metropolis" or algorithm_name == "hxy-metropolis" or
         algorithm_name == "xy-gaussian-noise-metropolis" or algorithm_name == "hxy-gaussian-noise-metropolis" or
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
         algorithm_name == "hxy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis" or
         algorithm_name == "hxy-gaussian-noise-metropolis") and (
            summary_statistic_string == "inverse_permittivity" or
            summary_statistic_string == "topological_sector_fluctuations")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity or "
              "topological_sector_fluctuations as the second positional argument.")
        exit()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
