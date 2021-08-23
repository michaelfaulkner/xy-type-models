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


def main(config_file, observable_string, no_of_trispectrum_octaves=4, trispectrum_base_period_shift=1):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    config_data_getter.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = config_data_getter.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)

    if no_of_jobs > 1:
        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)
    else:
        pool = None

    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"

        power_trispectrum = polyspectra.try_to_load_normalised_power_trispectrum(
            observable_string, output_directory, temperature_directory, beta, no_of_sites, no_of_equilibration_sweeps,
            no_of_jobs, pool, trispectrum_base_period_shift, no_of_trispectrum_octaves)

        figure, axis = plt.subplots(1, 2, figsize=(10, 5))
        [axis[index].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10) for index in range(2)]
        axis[0].set_ylabel(fr"$|S_X^3 \left( f, f' \right)|$ / $|S_X^3 \left( f_0, f' \right)|$", fontsize=10,
                           labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

        colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_trispectrum_octaves + 2)))
        for index in range(no_of_trispectrum_octaves + 2):
            current_color = next(colors)
            if index == 0:
                axis[0].loglog(power_trispectrum[1], power_trispectrum[2][index], color=current_color, label=r"f' = 0")
            else:
                axis[0].loglog(power_trispectrum[1], power_trispectrum[2][index], color=current_color,
                               label=fr"f' = {2 ** index}" " x " fr"{power_trispectrum[0][0]:.2e}")

        colors = iter(reversed(plt.cm.rainbow(np.linspace(0, 1, no_of_trispectrum_octaves + 2))))
        for index in range(no_of_trispectrum_octaves + 1, -1, -1):
            current_color = next(colors)
            if index == 0:
                axis[1].loglog(power_trispectrum[1], power_trispectrum[2][index], color=current_color, label=r"f' = 0")
            else:
                axis[1].loglog(power_trispectrum[1], power_trispectrum[2][index], color=current_color,
                               label=fr"f' = {2 ** index}" " x " fr"{power_trispectrum[0][0]:.2e}")

        figure.tight_layout()
        trispectrum_legend = (axis[0].legend(loc="lower left", fontsize=7.5),
                              axis[1].legend(loc="lower left", fontsize=7.5))
        for legend in trispectrum_legend:
            legend.get_frame().set_edgecolor("k")
            legend.get_frame().set_lw(1.5)
        figure.savefig(f"{output_directory}/{observable_string}_compare_power_trispectrum_frequencies_"
                       f"{no_of_trispectrum_octaves}_octaves_temp_eq_{temperature:.2f}.pdf", bbox_inches="tight")
        figure.clf()
        temperature -= magnitude_of_temperature_increments

    if no_of_jobs > 1:
        pool.close()


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate.  In addition, you "
                        "may provide no_of_trispectrum_octaves (default value is 4) and trispectrum_base_period_shift "
                        "(default value is 1) in the third and fourth positions (respectively).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  In addition, you may provide no_of_trispectrum_octaves (default "
              "value is 4) and trispectrum_base_period_shift (default value is 1) in the third and fourth positions "
              "(respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The third must be no_of_trispectrum_octaves.  In addition, you may"
              " provide trispectrum_base_period_shift (default value is 1) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The third / fourth must be no_of_trispectrum_octaves / "
              "trispectrum_base_period_shift.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
