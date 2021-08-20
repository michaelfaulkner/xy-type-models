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


def main(config_file, observable_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    check_for_config_errors(algorithm_name, no_of_jobs, observable_string)
    temperature = final_temperature

    beta = 1.0 / temperature
    temperature_directory = f"temp_eq_{temperature:.2f}"
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
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


def check_for_config_errors(algorithm_name, no_of_jobs, observable_string):
    if (observable_string != "helicity_modulus" and observable_string != "magnetisation_norm" and
            observable_string != "magnetisation_phase" and observable_string != "specific_heat" and
            observable_string != "inverse_permittivity" and observable_string != "topological_sector_fluctuations"):
        print("ConfigurationError: Give one of helicity_modulus, magnetisation_norm, specific_heat, "
              "inverse_permittivity or topological_sector_fluctuations as the second positional argument.")
        raise SystemExit
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (observable_string == "magnetisation_norm" or observable_string == "magnetisation_phase" or
             observable_string == "helicity_modulus")):
        print("ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase or helicity_modulus as the second positional argument.")
        raise SystemExit
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-metropolis" or
         algorithm_name == "hxy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis" or
         algorithm_name == "hxy-gaussian-noise-metropolis") and (
            observable_string == "inverse_permittivity" or observable_string == "topological_sector_fluctuations")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity or "
              "topological_sector_fluctuations as the second positional argument.")
        raise SystemExit
    if no_of_jobs != 1:
        print("ConfigurationError: Give a configuration file whose value of no_of_jobs is equal to one.")
        raise SystemExit


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("IndexError: Two positional arguments required - give the configuration-file location and "
              "the string of the observable whose trace plot you wish to create.")
        raise SystemExit
    main(sys.argv[1], sys.argv[2])
