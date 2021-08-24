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


def main(config_file, observable_string, no_of_power_2_correlators=3, no_of_power_10_correlators=4):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    config_data_getter.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = config_data_getter.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)

    figure, axis = plt.subplots(1 + no_of_power_2_correlators + no_of_power_10_correlators, figsize=(10, 20))
    plt.xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))

    if no_of_jobs > 1:
        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)
    else:
        pool = None

    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"

        power_spectrum = polyspectra.get_normalised_power_spectrum(observable_string, output_directory,
                                                                   temperature_directory, beta, no_of_sites,
                                                                   no_of_equilibration_sweeps, no_of_jobs, pool)
        power_spectrum_of_correlators = []
        for index in range(no_of_power_2_correlators):
            power_spectrum_of_correlators.append(polyspectra.get_power_spectrum_of_correlator(
                observable_string, output_directory, temperature_directory, beta, no_of_sites,
                no_of_equilibration_sweeps, no_of_jobs, pool, 2 ** (index + 1)))
        for index in range(no_of_power_10_correlators):
            power_spectrum_of_correlators.append(polyspectra.get_power_spectrum_of_correlator(
                observable_string, output_directory, temperature_directory, beta, no_of_sites,
                no_of_equilibration_sweeps, no_of_jobs, pool, 10 ** (index + 1)))

        current_color = next(colors)
        axis[0].loglog(power_spectrum[0], power_spectrum[1], color=current_color)
        for index, spectrum in enumerate(power_spectrum_of_correlators):
            if index == 0:
                axis[index + 1].loglog(spectrum[0], spectrum[1], color=current_color,
                                       label=fr"temperature = {temperature:.2f}, $\Delta t = "
                                             fr"{0.5 / power_spectrum[0, 0] / len(power_spectrum[0]):.2e}$")
            else:
                axis[index + 1].loglog(spectrum[0], spectrum[1], color=current_color)

        temperature -= magnitude_of_temperature_increments

    if no_of_jobs > 1:
        pool.close()

    x = np.linspace(1.0e-3, 1.0, 10000)
    axis[0].loglog(x, 1.0e-3 * x ** (-1.0), color="red", label=r"$f^{-1}$")
    axis[0].loglog(x, 5.0e-5 * x ** (-1.4), color="black", label=r"$f^{-1.4}$")
    axis[0].set_ylabel(r"$S_X \left( f \right)$ / $S_X \left( f_0 \right)$", fontsize=10, labelpad=10)
    for index in range(no_of_power_2_correlators + no_of_power_10_correlators):
        if index < no_of_power_2_correlators:
            axis[index + 1].set_ylabel(fr"$S_Y \left( f \right)$ / $S_Y \left( f_0 \right)$, $Y(t) = X(t) "
                                       fr"X(t + {2 ** (index + 1)} \Delta t)$", fontsize=7.5, labelpad=10)
        else:
            axis[index + 1].set_ylabel(fr"$S_Y \left( f \right)$ / $S_Y \left( f_0 \right)$, $Y(t) = X(t) "
                                       fr"X(t + {10 ** (index - no_of_power_2_correlators + 1)} \Delta t)$",
                                       fontsize=7.5, labelpad=10)
    figure.tight_layout()
    correlators_legend = [axis[index].legend(loc="lower left", fontsize=10) for index in range(2)]
    for legend in correlators_legend:
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/{observable_string}_power_spectra_of_signal_and_correlators.pdf",
                   bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate.  In addition, you "
                        "may provide no_of_power_2_correlators (default value is 3) and no_of_power_10_correlators "
                        "(default value is 4) in the third and fourth positions (respectively).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  In addition, you may provide no_of_power_2_correlators (default "
              "value is 3) and no_of_power_10_correlators (default value is 4) in the third and fourth positions "
              "(respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The third must be no_of_power_2_correlators.  In addition, you may"
              " provide no_of_power_10_correlators (default value is 4) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The third / fourth must be no_of_power_2_correlators / "
              "no_of_power_10_correlators.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
