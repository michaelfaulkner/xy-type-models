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
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    config_data_getter.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = config_data_getter.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)

    no_of_trispectrum_octaves = 6
    trispectrum_base_period_shift = 1

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
