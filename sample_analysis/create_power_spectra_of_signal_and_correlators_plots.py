import importlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import time

# import additional modules
setup_scripts = importlib.import_module("setup_scripts")
polyspectra = importlib.import_module("polyspectra")


def main(config_file, observable_string, no_of_power_2_correlators=3, no_of_power_10_correlators=4):
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
     no_of_jobs, temperature, magnitude_of_temperature_increments, pool) = setup_scripts.set_up_polyspectra_script(
        config_file, observable_string)

    figure, axes = plt.subplots(2 + no_of_power_2_correlators + no_of_power_10_correlators, figsize=(10, 20))
    plt.xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))

    start_time = time.time()
    for i in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")

        power_spectrum = polyspectra.get_power_spectrum(algorithm_name, observable_string, output_directory,
                                                        temperature, no_of_sites, no_of_equilibration_sweeps,
                                                        no_of_jobs, pool)
        # normalise power spectrum with respect to its low-frequency values
        power_spectrum[1] /= power_spectrum[1, 0]

        second_spectrum = polyspectra.get_second_spectrum(
            algorithm_name, observable_string, output_directory, temperature, no_of_sites, no_of_equilibration_sweeps,
            no_of_jobs, pool)
        # normalise second spectrum spectrum with respect to its low-frequency value
        second_spectrum[1] /= second_spectrum[1, 0]

        power_spectra_of_correlators = []
        for index in range(no_of_power_2_correlators):
            power_spectrum_of_correlator = polyspectra.get_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, no_of_jobs, pool, 2 ** index)
            # normalise power spectrum with respect to its low-frequency value
            power_spectrum_of_correlator[1] /= power_spectrum_of_correlator[1, 0]
            power_spectra_of_correlators.append(power_spectrum_of_correlator)
        for index in range(no_of_power_10_correlators):
            power_spectrum_of_correlator = polyspectra.get_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, no_of_sites,
                no_of_equilibration_sweeps, no_of_jobs, pool, 10 ** (index + 1))
            # normalise power spectrum with respect to its low-frequency value
            power_spectrum_of_correlator[1] /= power_spectrum_of_correlator[1, 0]
            power_spectra_of_correlators.append(power_spectrum_of_correlator)

        current_color = next(colors)
        max_power_spectrum_index = np.argmax(power_spectrum[0] > 1.0e-1) - 1
        axes[0].loglog(power_spectrum[0, :max_power_spectrum_index], power_spectrum[1, :max_power_spectrum_index],
                       color=current_color, label=f"temperature = {temperature:.2f}")
        max_second_spectrum_index = np.argmax(second_spectrum[0] > 1.0e-1) - 1
        axes[1].loglog(second_spectrum[0, :max_second_spectrum_index], second_spectrum[1, :max_second_spectrum_index],
                       color=current_color, label=f"temperature = {temperature:.2f}")
        for index, spectrum in enumerate(power_spectra_of_correlators):
            max_spectrum_index = np.argmax(spectrum[0] > 1.0e-1) - 1
            axes[index + 2].loglog(spectrum[0, :max_spectrum_index], spectrum[1, :max_spectrum_index],
                                   color=current_color,
                                   label=fr"temperature = {temperature:.2f}, $\Delta t = "
                                         fr"{0.5 / power_spectrum[0, 0] / len(power_spectrum[0]):.2e}$")

        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    if no_of_jobs > 1:
        pool.close()

    axes[0].set_ylabel(r"$S_X(f)$ / $S_X(f_0)$", fontsize=10, labelpad=10)
    axes[1].set_ylabel(r"$S_Y(f)$ / $S_Y(f_0)$, $Y(t) = X(t)X(t)$", fontsize=10, labelpad=10)
    for index in range(no_of_power_2_correlators + no_of_power_10_correlators):
        if index < no_of_power_2_correlators:
            axes[index + 2].set_ylabel(fr"$S_Y(f)$ / $S_Y(f_0)$, $Y(t) = X(t) "
                                       fr"X(t + {2 ** index} \Delta t)$", fontsize=7.5, labelpad=10)
        else:
            axes[index + 2].set_ylabel(fr"$S_Y(f)$ / $S_Y(f_0)$, $Y(t) = X(t) "
                                       fr"X(t + {10 ** (index - no_of_power_2_correlators + 1)} \Delta t)$",
                                       fontsize=7.5, labelpad=10)

    figure.tight_layout()
    correlators_legend = [axis.legend(loc="lower left", fontsize=10) for axis in axes]
    for legend in correlators_legend:
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/{observable_string}_power_spectra_of_signal_and_correlators_"
                   f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                   bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose polyspectra you wish to estimate (in the first and second "
                        "positions, respectively).  In addition, you may provide no_of_power_2_correlators (default "
                        "value is 3) and no_of_power_10_correlators (default value is 4) in the third and fourth "
                        "positions (respectively).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file / "
              "the string of the observable whose polyspectra you wish to estimate.  In addition, you may provide "
              "no_of_power_2_correlators (default value is 3) and no_of_power_10_correlators (default value is 4) in "
              "the third and fourth positions (respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file / the string of the observable whose polyspectra you wish to estimate / "
              "no_of_power_2_correlators.  In addition, you may provide no_of_power_10_correlators (default value is 4)"
              " in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The first / second / third / fourth must be the location of the "
              "configuration file / the string of the observable whose polyspectra you wish to estimate / "
              "no_of_power_2_correlators / no_of_power_10_correlators.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
