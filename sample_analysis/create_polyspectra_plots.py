from scipy.optimize import curve_fit
import importlib
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
polyspectra = importlib.import_module("polyspectra")
setup_scripts = importlib.import_module("setup_scripts")


def main(config_file, observable_string, no_of_trispectrum_auxiliary_frequency_octaves=6,
         trispectrum_base_period_shift=1):
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
     no_of_jobs, temperature, magnitude_of_temperature_increments, pool) = setup_scripts.set_up_polyspectra_script(
        config_file, observable_string)

    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))
    figure, axes = plt.subplots(2, figsize=(10, 10))
    axes[1].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    axes[0].set_ylabel(r"$S_X \left( f \right)$ / $S_X \left( f_0 \right)$", fontsize=10, labelpad=10)
    axes[1].set_ylabel(r"$\left| S_X^3 \left( f, f' \right) \right|$ / $\left| S_X^3 \left( f_0, f' \right) \right|$",
                       fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

    start_time = time.time()
    for temperature_index in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")
        current_color = next(colors)

        power_spectrum = polyspectra.get_power_spectrum(algorithm_name, observable_string, output_directory,
                                                        temperature, no_of_sites, no_of_equilibration_sweeps,
                                                        no_of_jobs, pool)
        power_trispectrum = polyspectra.get_power_trispectrum(
            algorithm_name, observable_string, output_directory, temperature, no_of_sites, no_of_equilibration_sweeps,
            no_of_jobs, pool, no_of_trispectrum_auxiliary_frequency_octaves, trispectrum_base_period_shift)

        # normalise polyspectra with respect to their low-frequency values
        power_spectrum[1] /= power_spectrum[1, 0]
        power_trispectrum[2] = [spectrum / spectrum[0] for spectrum in power_trispectrum[2]]

        # note frequency indices in order to discard spectrum (trispectrum) for all frequencies >= 0.1 (5.0e-3)  since
        # i) interesting physics at long timescales, and ii) emergent Langevin diffusion breaks down at short timescales
        max_power_spectrum_index = np.argmax(power_spectrum[0] > 1.0e-1) - 1
        max_power_trispectrum_index = np.argmax(power_trispectrum[1] > 5.0e-3) - 1

        """fit Lorentzian model to power spectrum"""
        initial_frequency = power_spectrum[0, 0]
        final_frequency = 1.0e-2
        increment = 10.0 ** math.floor(math.log(initial_frequency, 10))
        final_frequency_index = np.argmax(power_spectrum[0] > final_frequency) - 1
        if no_of_jobs > 1:
            parameter_values_and_errors = curve_fit(lorentzian_model, power_spectrum[0, :final_frequency_index],
                                                    power_spectrum[1, :final_frequency_index],
                                                    sigma=power_spectrum[2, :final_frequency_index],
                                                    bounds=(np.array([0.1, initial_frequency]),
                                                            np.array([10.0, final_frequency])))
        else:
            parameter_values_and_errors = curve_fit(lorentzian_model, power_spectrum[0, :final_frequency_index],
                                                    power_spectrum[1, :final_frequency_index],
                                                    bounds=(np.array([0.1, initial_frequency]),
                                                            np.array([10.0, final_frequency])))
        parameter_values = parameter_values_and_errors[0]
        parameter_errors = np.sqrt(np.diag(parameter_values_and_errors[1]))
        model_frequency_values = np.arange(initial_frequency, final_frequency, increment)
        model_spectrum_values = lorentzian_model(model_frequency_values, *parameter_values)
        axes[0].loglog(power_spectrum[0, :max_power_spectrum_index], power_spectrum[1, :max_power_spectrum_index],
                       color=current_color,
                       label=fr"temperature = {temperature:.2f}; $f_c$ = {parameter_values[1]:.2e} $\pm$ "
                             fr"{parameter_errors[1]:.2e}; $S_0$ = {parameter_values[0]:.2f} $\pm$ "
                             fr"{parameter_errors[0]:.2e}")
        axes[0].loglog(model_frequency_values, model_spectrum_values, color='k')

        """lower frequency range of trispectrum fitting"""
        (one_over_f_model_lower_parameters, one_over_f_model_lower_errors, one_over_f_model_lower_frequency_values,
         one_over_f_model_lower_spectrum_values) = fit_one_over_f_model_to_trispectrum(power_trispectrum, "lower", 10.0,
                                                                                       no_of_sites, no_of_jobs)
        """upper frequency range of trispectrum fitting"""
        (one_over_f_model_upper_parameters, one_over_f_model_upper_errors, one_over_f_model_upper_frequency_values,
         one_over_f_model_upper_spectrum_values) = fit_one_over_f_model_to_trispectrum(power_trispectrum, "upper", 10.0,
                                                                                       no_of_sites, no_of_jobs)
        axes[1].loglog(power_trispectrum[1][:max_power_trispectrum_index],
                       power_trispectrum[2][len(power_trispectrum[2]) - 1][:max_power_trispectrum_index],
                       color=current_color,
                       label=fr"temperature = {temperature:.2f}; "
                             fr"f' = {power_trispectrum[0][len(power_trispectrum[0]) - 1]:.2e}; "
                             fr"$\alpha_L$ = {one_over_f_model_lower_parameters[1]:.2e} $\pm$ "
                             fr"{one_over_f_model_lower_errors[1]:.2e}; "
                             fr"$\alpha_U$ = {one_over_f_model_upper_parameters[1]:.2e} $\pm$ "
                             fr"{one_over_f_model_upper_errors[1]:.2e}")
        axes[1].loglog(one_over_f_model_lower_frequency_values, one_over_f_model_lower_spectrum_values, color='k')
        axes[1].loglog(one_over_f_model_upper_frequency_values, one_over_f_model_upper_spectrum_values, color='k',
                       linestyle='dashed')

        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    figure.tight_layout()
    if observable_string == "cartesian_magnetisation":
        legends = [axes[0].legend(title=r"black lines: fits to $S_X(f)$ / $S_X(f_0) = S_0 / (1 + (f / f_c)^2)$",
                                  loc="lower left", fontsize=10),
                   axes[1].legend(title=r"black lines: fits to $S_X^3(f, f') \sim f^{-\alpha}$; $\alpha = \alpha_L$ "
                                        r"(solid), $\alpha = \alpha_U$ (dashed)", loc="lower left", fontsize=10)]
    else:
        legends = [axis.legend(loc="lower left", fontsize=10) for axis in axes]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(1.5) for legend in legends]
    figure.savefig(f"{output_directory}/{observable_string}_polyspectra_max_trispectrum_shift_eq_"
                   f"{2 ** no_of_trispectrum_auxiliary_frequency_octaves}_x_{trispectrum_base_period_shift}_delta_t_"
                   f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                   bbox_inches="tight")
    if no_of_jobs > 1:
        pool.close()


def fit_one_over_f_model_to_trispectrum(power_trispectrum, frequency_range, max_model_exponent, no_of_sites,
                                        no_of_jobs):
    """fit one-over-f model to power trispectrum"""
    if frequency_range == "lower":
        if no_of_sites >= 64 ** 2:
            initial_frequency, final_frequency = 4.0e-5, 1.5e-4
        elif no_of_sites == 32 ** 2:
            initial_frequency, final_frequency = 1.2e-4, 2.0e-4
        elif no_of_sites <= 16 ** 2:
            initial_frequency, final_frequency = 4.0e-4, 8.0e-4
    elif frequency_range == "upper":
        if no_of_sites >= 64 ** 2:
            initial_frequency, final_frequency = 1.5e-4, 4.0e-4
        elif no_of_sites == 32 ** 2:
            initial_frequency, final_frequency = 2.0e-4, 6.0e-4
        elif no_of_sites <= 16 ** 2:
            initial_frequency, final_frequency = 8.0e-4, 1.5e-3
    else:
        raise SystemExit("Give 'lower' or 'upper' as the value for frequency_range in the "
                         "fit_one_over_f_model_to_trispectrum() method.")
    increment = 10.0 ** math.floor(math.log(initial_frequency, 10))
    initial_frequency_index = np.argmax(power_trispectrum[1] > initial_frequency) - 1
    final_frequency_index = np.argmax(power_trispectrum[1] > final_frequency) - 1
    # the following commented-out code may be useful for no_of_jobs very large; we found worse results using the
    # trispectrum error bars for no_of_jobs = 8
    """if no_of_jobs > 1:
        parameter_values_and_errors = curve_fit(
            one_over_f_model, power_trispectrum[1][initial_frequency_index:final_frequency_index],
            power_trispectrum[2][len(power_trispectrum[2]) - 1][initial_frequency_index:final_frequency_index],
            sigma=power_trispectrum[3][len(power_trispectrum[3]) - 1][initial_frequency_index:final_frequency_index],
            bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))
    else:
        parameter_values_and_errors = curve_fit(
            one_over_f_model, power_trispectrum[1][initial_frequency_index:final_frequency_index],
            power_trispectrum[2][len(power_trispectrum[2]) - 1][initial_frequency_index:final_frequency_index],
            bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))"""
    parameter_values_and_errors = curve_fit(
        one_over_f_model, power_trispectrum[1][initial_frequency_index:final_frequency_index],
        power_trispectrum[2][len(power_trispectrum[2]) - 1][initial_frequency_index:final_frequency_index],
        bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))
    parameter_values = parameter_values_and_errors[0]
    parameter_errors = np.sqrt(np.diag(parameter_values_and_errors[1]))
    model_frequency_values = np.arange(initial_frequency, final_frequency, increment)
    model_spectrum_values = one_over_f_model(model_frequency_values, *parameter_values)
    return parameter_values, parameter_errors, model_frequency_values, model_spectrum_values


def lorentzian_model(frequencies, zero_frequency_value, characteristic_frequency):
    return zero_frequency_value / (1.0 + (frequencies / characteristic_frequency) ** 2.0)


def one_over_f_model(frequencies, scale_factor, exponent):
    return scale_factor / frequencies ** exponent


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate (in the first and "
                        "second positions, respectively).  In addition, you may provide "
                        "no_of_trispectrum_auxiliary_frequency_octaves (default value is 6) and "
                        "trispectrum_base_period_shift (default value is 1) in the third and fourth positions, "
                        "respectively.")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file / "
              "the string of the observable whose power trispectrum you wish to estimate.  In addition, you may provide"
              " no_of_trispectrum_auxiliary_frequency_octaves (default value is 6) and trispectrum_base_period_shift "
              "(default value is 1) in the third and fourth positions (respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file / the string of the observable whose power trispectrum you wish to estimate / "
              "no_of_trispectrum_auxiliary_frequency_octaves.  In addition, you may provide "
              "trispectrum_base_period_shift (default value is 1) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The first / second / third / fourth must be the location of the "
              "configuration file / the string of the observable whose power trispectrum you wish to estimate / "
              "no_of_trispectrum_auxiliary_frequency_octaves / trispectrum_base_period_shift.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
