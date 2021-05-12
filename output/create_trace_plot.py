import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")


def main(config_file, summary_statistic_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    basic_config_data = config_data_getter.get_basic_data(config_file)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, temperature, no_of_jobs) = (
        basic_config_data[0], basic_config_data[1], basic_config_data[2], basic_config_data[3], basic_config_data[5],
        basic_config_data[8])

    check_for_config_errors(algorithm_name, no_of_jobs, summary_statistic_string)
    no_of_sites = integer_lattice_length ** 2
    beta = 1.0 / temperature
    temperature_directory = f"temp_eq_{temperature:.2f}"
    get_sample_method = getattr(sample_getter, "get_" + summary_statistic_string)
    sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)

    plt.xlabel(r"time, $t$", fontsize=15, labelpad=10)
    plt.ylabel(r"$ x \left( t \right)$", fontsize=15, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
    plt.plot(sample, color="k", linewidth=1, linestyle="-")
    plt.tight_layout()
    plt.show()
    plt.clf()
    plt.plot(sample[0:no_of_equilibration_sweeps], color="k", linewidth=1, linestyle="-")
    plt.tight_layout()
    plt.show()


def check_for_config_errors(algorithm_name, no_of_jobs, summary_statistic_string):
    if (summary_statistic_string != "helicity_modulus" and summary_statistic_string != "magnetisation_norm" and
            summary_statistic_string != "magnetisation_phase" and summary_statistic_string != "specific_heat" and
            summary_statistic_string != "inverse_permittivity" and
            summary_statistic_string != "topological_sector_fluctuations"):
        print("ConfigurationError: Give one of helicity_modulus, magnetisation_norm, specific_heat, "
              "inverse_permittivity or topological_sector_fluctuations as the second positional argument.")
        exit()
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (summary_statistic_string == "magnetisation_norm" or summary_statistic_string == "magnetisation_phase" or
             summary_statistic_string == "helicity_modulus")):
        print("ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase or helicity_modulus as the second positional argument.")
        exit()
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-metropolis" or
         algorithm_name == "hxy-metropolis") and (summary_statistic_string == "inverse_permittivity" or
                                                  summary_statistic_string == "topological_sector_fluctuations")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity or "
              "topological_sector_fluctuations as the second positional argument.")
        exit()
    if no_of_jobs != 1:
        print("ConfigurationError: Give a configuration file whose value of no_of_jobs is equal to one.")
        exit()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
