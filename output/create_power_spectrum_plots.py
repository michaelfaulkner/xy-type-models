import importlib
import matplotlib
import matplotlib.pyplot as plt
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
polyspectra = importlib.import_module("polyspectra")


def main(config_file, power_spectrum_string):
    basic_config_data = config_data_getter.get_basic_data(config_file)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = (basic_config_data[0], basic_config_data[1],
                                                                     basic_config_data[2], basic_config_data[3],
                                                                     basic_config_data[5], basic_config_data[6],
                                                                     basic_config_data[7], basic_config_data[8])

    check_for_config_errors(algorithm_name, power_spectrum_string)
    no_of_sites = integer_lattice_length ** 2
    temperature = final_temperature
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / no_of_temperature_increments

    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    figure, axis = plt.subplots(2)
    plt.xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    plt.ylabel(r"$ S_X \left( f \right)$", fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))

    if no_of_jobs > 1:
        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)

    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"

        try:
            with open(f"{output_directory}/{power_spectrum_string}_power_spectrum_temp_eq_{temperature:.2f}.csv",
                      "r") as power_spectrum_file:
                power_spectrum = np.loadtxt(power_spectrum_file, dtype=float, delimiter=",")
        except IOError:
            if no_of_jobs == 1:
                power_spectrum = polyspectra.get_power_spectrum(power_spectrum_string, output_directory,
                                                                temperature_directory, beta, no_of_sites,
                                                                no_of_equilibration_sweeps)
            else:
                power_spectra = pool.starmap(polyspectra.get_power_spectrum,
                                             [(power_spectrum_string, f"{output_directory}/job_{job_number + 1}",
                                               temperature_directory, beta, no_of_sites, no_of_equilibration_sweeps)
                                              for job_number in range(no_of_jobs)])
                power_spectrum = np.mean(np.array(power_spectra), axis=0)
            with open(f"{output_directory}/{power_spectrum_string}_power_spectrum_temp_eq_{temperature:.2f}.csv",
                      "w") as power_spectrum_file:
                np.savetxt(power_spectrum_file, power_spectrum, delimiter=",")
        current_color = next(colors)
        axis[0].plot(power_spectrum[0, 1:], power_spectrum[1, 1:], color=current_color,
                     label=f"temperature = {temperature:.2f}")
        axis[1].loglog(power_spectrum[0, 1:], power_spectrum[1, 1:], color=current_color)
        temperature -= magnitude_of_temperature_increments

    pool.close()
    x = np.linspace(1.0e-3, 10.0, 10000)
    y = 0.5 * x ** (-1.0)
    axis[1].loglog(x, y, color="black")
    y = 0.01 * x ** (-1.4)
    axis[1].loglog(x, y, color="black")
    axis[0].set_xlim(0.0, 0.5)
    axis[0].set_ylim(0.0, 5.0)
    axis[1].set_xlim(1.0e-5, 1000.0)
    axis[1].set_ylim(1.0e-5, 1000.0)
    figure.tight_layout()
    legend = axis[0].legend(loc="upper right", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/{power_spectrum_string}_power_spectrum.pdf", bbox_inches="tight")


def check_for_config_errors(algorithm_name, power_spectrum_string):
    if ((algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            (power_spectrum_string == "magnetisation_norm" or power_spectrum_string == "magnetisation_norm" or
             power_spectrum_string == "helicity_modulus" or power_spectrum_string == "inverse_vacuum_permittivity" or
             power_spectrum_string == "toroidal_vortex_polarisation")):
        print("ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation_norm, "
              "magnetisation_phase, helicity_modulus, inverse_vacuum_permittivity or toroidal_vortex_polarisation as "
              "the second positional argument.")
        exit()
    if ((algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc" or algorithm_name == "xy-metropolis" or
         algorithm_name == "hxy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis" or
         algorithm_name == "hxy-gaussian-noise-metropolis") and (
            power_spectrum_string == "inverse_permittivity" or
            power_spectrum_string == "topological_sector_fluctuations" or
            power_spectrum_string == "toroidal_polarisation")):
        print("ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity, "
              "topological_sector_fluctuations or toroidal_polarisation as the second positional argument.")
        exit()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("IndexError: Two positional arguments required - give the configuration-file location and "
              "the summary-statistic string whose power spectrum you wish to calculate.")
    main(sys.argv[1], sys.argv[2])
