from scipy.optimize import curve_fit
import importlib
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
import time

# import additional modules
setup_scripts = importlib.import_module("setup_scripts")
polyspectra = importlib.import_module("polyspectra")


def main(config_file, observable_string, no_of_power_2_correlators=3, no_of_power_10_correlators=4):
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
     external_global_moves_string, no_of_jobs, temperature, magnitude_of_temperature_increments,
     pool) = setup_scripts.set_up_polyspectra_script(config_file, observable_string)

    figure, axes = plt.subplots(2 + no_of_power_2_correlators + no_of_power_10_correlators, figsize=(10, 20))
    plt.xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))

    start_time = time.time()
    for temperature_index in reversed(range(no_of_temperature_increments + 1)):
        print(f"Temperature = {temperature:.2f}")

        power_spectrum = polyspectra.get_power_spectrum(algorithm_name, observable_string, output_directory,
                                                        temperature, temperature_index, no_of_sites,
                                                        no_of_equilibration_sweeps, external_global_moves_string,
                                                        no_of_jobs, pool)
        # normalise power spectrum with respect to its low-frequency values
        power_spectrum[1] /= power_spectrum[1, 0]

        second_spectrum = polyspectra.get_second_spectrum(
            algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
            no_of_equilibration_sweeps, external_global_moves_string, no_of_jobs, pool)
        # normalise second spectrum spectrum with respect to its low-frequency value
        second_spectrum[1] /= second_spectrum[1, 0]

        power_spectra_of_correlators = []
        for index in range(no_of_power_2_correlators):
            power_spectrum_of_correlator = polyspectra.get_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                no_of_equilibration_sweeps, external_global_moves_string, no_of_jobs, pool, 2 ** index)
            # normalise power spectrum with respect to its low-frequency value
            power_spectrum_of_correlator[1] /= power_spectrum_of_correlator[1, 0]
            power_spectra_of_correlators.append(power_spectrum_of_correlator)
        for index in range(no_of_power_10_correlators):
            power_spectrum_of_correlator = polyspectra.get_power_spectrum_of_correlator(
                algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
                no_of_equilibration_sweeps, external_global_moves_string, no_of_jobs, pool, 10 ** (index + 1))
            # normalise power spectrum with respect to its low-frequency value
            power_spectrum_of_correlator[1] /= power_spectrum_of_correlator[1, 0]
            power_spectra_of_correlators.append(power_spectrum_of_correlator)

        current_color = next(colors)

        (one_over_f_model_parameters, one_over_f_model_errors, one_over_f_model_frequency_values,
         one_over_f_model_spectrum_values) = fit_one_over_f_model_to_spectrum(power_spectrum, 10.0, no_of_jobs)
        max_power_spectrum_index = np.argmax(power_spectrum[0] > 1.0e-1) - 1
        axes[0].loglog(power_spectrum[0, :max_power_spectrum_index], power_spectrum[1, :max_power_spectrum_index],
                       color=current_color,
                       label=fr"temperature = {temperature:.2f}; "fr"$\alpha$ = {one_over_f_model_parameters[1]:.2e} "
                             fr"$\pm$ {one_over_f_model_errors[1]:.2e}")
        axes[0].loglog(one_over_f_model_frequency_values, one_over_f_model_spectrum_values, color='k')

        (one_over_f_model_parameters, one_over_f_model_errors, one_over_f_model_frequency_values,
         one_over_f_model_spectrum_values) = fit_one_over_f_model_to_spectrum(second_spectrum, 10.0, no_of_jobs)
        max_second_spectrum_index = np.argmax(second_spectrum[0] > 1.0e-1) - 1
        axes[1].loglog(second_spectrum[0, :max_second_spectrum_index], second_spectrum[1, :max_second_spectrum_index],
                       color=current_color,
                       label=fr"temperature = {temperature:.2f}; "fr"$\alpha$ = {one_over_f_model_parameters[1]:.2e} "
                             fr"$\pm$ {one_over_f_model_errors[1]:.2e}")
        axes[1].loglog(one_over_f_model_frequency_values, one_over_f_model_spectrum_values, color='k')

        for index, spectrum in enumerate(power_spectra_of_correlators):
            (one_over_f_model_parameters, one_over_f_model_errors, one_over_f_model_frequency_values,
             one_over_f_model_spectrum_values) = fit_one_over_f_model_to_spectrum(spectrum, 10.0, no_of_jobs)
            max_spectrum_index = np.argmax(spectrum[0] > 1.0e-1) - 1
            axes[index + 2].loglog(spectrum[0, :max_spectrum_index], spectrum[1, :max_spectrum_index],
                                   color=current_color,
                                   label=fr"temperature = {temperature:.2f}, $\Delta t = "
                                         fr"{0.5 / power_spectrum[0, 0] / len(power_spectrum[0]):.2e}$; "
                                         fr"$\alpha$ = {one_over_f_model_parameters[1]:.2e} $\pm$ "
                                         fr"{one_over_f_model_errors[1]:.2e}")
            axes[index + 2].loglog(one_over_f_model_frequency_values, one_over_f_model_spectrum_values, color='k')

        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    if no_of_jobs > 1:
        pool.close()

    axes[0].set_ylabel(r"$S_X(f)$ / $S_X(f_0)$", fontsize=10, labelpad=10)
    axes[1].set_ylabel(r"$S_Y(f)$ / $S_Y(f_0)$, $Y(t) = X(t)X(t)$", fontsize=10, labelpad=10)
    for index in range(no_of_power_2_correlators + no_of_power_10_correlators):
        if index < no_of_power_2_correlators:
            axes[index + 2].set_ylabel(fr"$S_Y(f)$ / $S_Y(f_0)$, $Y(t) = X(t) "
                                       fr"X(t - {2 ** index} \Delta t)$", fontsize=7.5, labelpad=10)
        else:
            axes[index + 2].set_ylabel(fr"$S_Y(f)$ / $S_Y(f_0)$, $Y(t) = X(t) "
                                       fr"X(t - {10 ** (index - no_of_power_2_correlators + 1)} \Delta t)$",
                                       fontsize=7.5, labelpad=10)

    figure.tight_layout()
    correlators_legend = [axis.legend(loc="lower left", fontsize=7.5) for axis in axes]
    for legend in correlators_legend:
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/{observable_string}_power_spectra_of_signal_and_correlators_"
                   f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                   f"{int(no_of_sites ** 0.5)}.pdf", bbox_inches="tight")


def fit_one_over_f_model_to_spectrum(spectrum, max_model_exponent, no_of_jobs):
    """fit one-over-f model to power spectrum"""
    initial_frequency, final_frequency = 5.0e-4, 3.0e-3
    increment = 10.0 ** math.floor(math.log(initial_frequency, 10)) / 2.0
    initial_frequency_index = np.argmax(spectrum[0] > initial_frequency) - 1
    final_frequency_index = np.argmax(spectrum[0] > final_frequency) - 1
    if no_of_jobs > 32:
        """n.b., no_of_jobs > 32 as we found worse results using the trispectrum error bars for no_of_jobs = 8."""
        parameter_values_and_errors = curve_fit(
            one_over_f_model, spectrum[0, initial_frequency_index:final_frequency_index],
            spectrum[1, initial_frequency_index:final_frequency_index],
            sigma=spectrum[2, initial_frequency_index:final_frequency_index],
            bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))
    else:
        parameter_values_and_errors = curve_fit(
            one_over_f_model, spectrum[0, initial_frequency_index:final_frequency_index],
            spectrum[1, initial_frequency_index:final_frequency_index],
            bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))
    parameter_values = parameter_values_and_errors[0]
    parameter_errors = np.sqrt(np.diag(parameter_values_and_errors[1]))
    model_frequency_values = np.arange(initial_frequency, final_frequency, increment)
    model_spectrum_values = one_over_f_model(model_frequency_values, *parameter_values)
    return parameter_values, parameter_errors, model_frequency_values, model_spectrum_values


def one_over_f_model(frequencies, scale_factor, exponent):
    return scale_factor / frequencies ** exponent


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
