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


def main(config_file, observable_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    config_data_getter.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = config_data_getter.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)

    no_of_power_2_correlators = 3
    no_of_power_10_correlators = 4
    no_of_trispectrum_octaves = 2
    trispectrum_base_period_shift = 1

    correlators_figure, correlators_axis = plt.subplots(1 + no_of_power_2_correlators + no_of_power_10_correlators,
                                                        figsize=(10, 20))
    plt.xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)
    trispectrum_figure, trispectrum_axis = plt.subplots(no_of_trispectrum_octaves + 2, figsize=(10, 10))
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
            compute_power_spectra_of_correlators(beta, index, 2, no_of_equilibration_sweeps, no_of_jobs, no_of_sites,
                                                 output_directory, pool, power_spectrum_of_correlators,
                                                 observable_string, temperature, temperature_directory)

        for index in range(no_of_power_10_correlators):
            compute_power_spectra_of_correlators(beta, index, 10, no_of_equilibration_sweeps, no_of_jobs, no_of_sites,
                                                 output_directory, pool, power_spectrum_of_correlators,
                                                 observable_string, temperature, temperature_directory)
        power_trispectrum = polyspectra.get_normalised_power_trispectrum(observable_string, output_directory,
                                                                         temperature_directory, beta, no_of_sites,
                                                                         no_of_equilibration_sweeps, no_of_jobs, pool,
                                                                         no_of_trispectrum_octaves,
                                                                         trispectrum_base_period_shift)

        current_color = next(colors)
        correlators_axis[0].loglog(power_spectrum[0], power_spectrum[1], color=current_color)
        for index, correlator in enumerate(power_spectrum_of_correlators):
            if index == 0:
                correlators_axis[index + 1].loglog(correlator[0], correlator[1], color=current_color,
                                                   label=f"temperature = {temperature:.2f}")
            else:
                correlators_axis[index + 1].loglog(correlator[0], correlator[1], color=current_color)
        for index in range(no_of_trispectrum_octaves + 2):
            if index == 0:
                trispectrum_axis[index].loglog(power_trispectrum[1], power_trispectrum[2][index], color=current_color,
                                               label=f"temperature = {temperature:.2f}")
            else:
                trispectrum_axis[index].loglog(power_trispectrum[1], power_trispectrum[2][index], color=current_color)

        temperature -= magnitude_of_temperature_increments

    if no_of_jobs > 1:
        pool.close()

    x = np.linspace(1.0e-2, 10.0, 10000)
    correlators_axis[0].loglog(x, 0.01 * x ** (-1.0), color="red", label=r"$f^{-1}$")
    correlators_axis[0].loglog(x, 0.001 * x ** (-1.4), color="black", label=r"$f^{-1.4}$")
    correlators_axis[0].set_ylabel(r"$S_X \left( f \right)$ / $S_X \left( f_0 \right)$", fontsize=10, labelpad=10)
    for index in range(no_of_power_2_correlators + no_of_power_10_correlators):
        if index < no_of_power_2_correlators:
            correlators_axis[index + 1].set_ylabel(
                fr"$S_Y \left( f \right)$ / $S_Y \left( f_0 \right)$, $Y(t) = X(t) X(t + {2 ** (index + 1)})$",
                fontsize=10,
                labelpad=10)
        else:
            correlators_axis[index + 1].set_ylabel(
                fr"$S_Y \left( f \right)$ / $S_Y \left( f_0 \right)$, "
                fr"$Y(t) = X(t) X(t + {10 ** (index - no_of_power_2_correlators + 1)})$", fontsize=10, labelpad=10)
    correlators_figure.tight_layout()
    correlators_legend = (correlators_axis[0].legend(loc="lower left", fontsize=10),
                          correlators_axis[1].legend(loc="lower left", fontsize=10))
    for legend in correlators_legend:
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(1.5)
    correlators_figure.savefig(f"{output_directory}/{observable_string}_normalised_power_spectrum.pdf",
                               bbox_inches="tight")

    for index in range(no_of_trispectrum_octaves + 2):
        if index == 0:
            trispectrum_axis[index].set_ylabel(fr"$S_X^3 \left( f, f' = 0 \right)$ / $S_X^3 \left( f_0, f' = 0 "
                                               fr"\right)$", fontsize=10, labelpad=10)
        elif index == 1:
            trispectrum_axis[index].set_ylabel(fr"$|S_X^3 \left( f, f_0' \right)|$ / $|S_X^3 \left( f_0, f_0' "
                                               fr"\right)|$, $f_0' = {power_trispectrum[0][0]:.2f}$", fontsize=10,
                                               labelpad=10)
        else:
            trispectrum_axis[index].set_ylabel(fr"$|S_X^3 \left( f, {2 ** (index - 1)} f_0' \right)|$ / $|S_X^3 \left( "
                                               fr"f_0, {2 ** (index - 1)} f_0' \right)|$, "
                                               fr"$f_0' = {power_trispectrum[0][0]:.2f}$", fontsize=10, labelpad=10)
    trispectrum_figure.tight_layout()
    trispectrum_legend = (correlators_axis[0].legend(loc="lower left", fontsize=10),
                          correlators_axis[1].legend(loc="lower left", fontsize=10))
    for legend in trispectrum_legend:
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(1.5)
    trispectrum_figure.savefig(f"{output_directory}/{observable_string}_normalised_power_trispectrum.pdf",
                               bbox_inches="tight")


def compute_power_spectra_of_correlators(beta, index, base, no_of_equilibration_sweeps, no_of_jobs, no_of_sites,
                                         output_directory, pool, power_spectrum_of_correlators, power_spectrum_string,
                                         temperature, temperature_directory):
    try:
        with open(f"{output_directory}/{power_spectrum_string}_normalised_power_spectrum_of_correlator_"
                  f"time_shift_eq_{base ** (index + 1)}_temp_eq_{temperature:.2f}.csv", "r") as data_file:
            correlator_power_spectrum = np.loadtxt(data_file, dtype=float, delimiter=",")
    except IOError:
        if no_of_jobs == 1:
            correlator_power_spectrum = polyspectra.get_power_spectrum_of_correlator(power_spectrum_string,
                                                                                     output_directory,
                                                                                     temperature_directory,
                                                                                     beta, no_of_sites,
                                                                                     no_of_equilibration_sweeps,
                                                                                     time_period_shift=(
                                                                                             base ** (index + 1)))
        else:
            correlator_power_spectra = pool.starmap(
                polyspectra.get_power_spectrum_of_correlator,
                [(power_spectrum_string, f"{output_directory}/job_{job_number + 1}", temperature_directory,
                  beta, no_of_sites, no_of_equilibration_sweeps, base ** (index + 1))
                 for job_number in range(no_of_jobs)])
            correlator_power_spectrum = np.mean(np.array(correlator_power_spectra), axis=0)
        # normalise power spectrum with respect to its low-frequency value
        correlator_power_spectrum[1] /= correlator_power_spectrum[1, 0]
        with open(f"{output_directory}/{power_spectrum_string}_normalised_power_spectrum_of_correlator_"
                  f"time_shift_eq_{base ** (index + 1)}_temp_eq_{temperature:.2f}.csv", "w") as data_file:
            np.savetxt(data_file, correlator_power_spectrum, delimiter=",")
    power_spectrum_of_correlators.append(correlator_power_spectrum)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose polyspectra you wish to estimate.")
    main(sys.argv[1], sys.argv[2])
