import importlib
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
polyspectra = importlib.import_module("polyspectra")


def main(config_file, observable_string):
    basic_config_data = config_data_getter.get_basic_data(config_file)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = (basic_config_data[0], basic_config_data[1],
                                                                     basic_config_data[2], basic_config_data[3],
                                                                     basic_config_data[5], basic_config_data[6],
                                                                     basic_config_data[7], basic_config_data[8])
    config_data_getter.check_for_config_errors(algorithm_name, observable_string)
    no_of_sites = integer_lattice_length ** 2
    temperature = final_temperature
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / no_of_temperature_increments
    no_of_trispectrum_octaves = 6
    trispectrum_base_period_shift = 1

    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    if no_of_jobs > 1:
        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)
    else:
        pool = None

    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"
        power_trispectrum = polyspectra.get_power_trispectrum(observable_string, output_directory,
                                                              temperature_directory, beta, no_of_sites,
                                                              no_of_equilibration_sweeps, no_of_jobs, pool,
                                                              no_of_trispectrum_octaves, trispectrum_base_period_shift)

        plt.xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
        plt.ylabel(fr"$|S_X^3 \left( f, f' \right)|$ / $|S_X^3 \left( f_0, f' \right)|$", fontsize=10, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=10, pad=10)
        plt.loglog(power_trispectrum[1], power_trispectrum[2][1], color='blue',
                   label=fr"f' = {power_trispectrum[0][0]}")
        plt.loglog(power_trispectrum[1], power_trispectrum[2][2], color='green',
                   label=r"f' = 2" " x " fr"{power_trispectrum[0][0]}")
        plt.loglog(power_trispectrum[1], power_trispectrum[2][len(power_trispectrum[2]) - 1], color='black',
                   label=fr"f' = {2 ** no_of_trispectrum_octaves}" " x " fr"{power_trispectrum[0][0]}")
        plt.loglog(power_trispectrum[1], power_trispectrum[2][0], color='red', label=r"f' = 0")
        plt.tight_layout()
        legend = plt.legend(loc="lower left", fontsize=10)
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(1.5)
        plt.savefig(f"{output_directory}/{observable_string}_compare_power_trispectrum_frequencies_"
                    f"{no_of_trispectrum_octaves}_octaves_temp_eq_{temperature:.2f}.pdf", bbox_inches="tight")
        plt.clf()

        temperature -= magnitude_of_temperature_increments

    if no_of_jobs > 1:
        pool.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("IndexError: Two positional arguments required - give the configuration-file location and "
              "the string of the observable whose power spectrum you wish to calculate.")
        raise SystemExit
    main(sys.argv[1], sys.argv[2])
