from polyspectra import get_power_spectrum, get_second_spectrum, get_power_trispectrum_nonzero_mode
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from setup_scripts import reverse_enumerate, setup_polyspectra_script
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import warnings


def main(config_file, observable_string, no_of_trispectrum_auxiliary_frequency_octaves=6,
         trispectrum_base_period_shift=1, target_auxiliary_frequency=None):
    (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, temperatures,
     external_global_moves_string, no_of_runs, pool) = setup_polyspectra_script(config_file, observable_string)

    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(temperatures) + 1)))
    figure, axes = plt.subplots(3, figsize=(10, 10))
    axes[2].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    axes[0].set_ylabel(r"$S_X \left( f \right)$ / $S_X \left( f_0 \right)$", fontsize=10, labelpad=10)
    axes[1].set_ylabel(r"$S_X^2 \left( f \right)$ / $S_X^2 \left( f_0 \right)$", fontsize=10, labelpad=10)
    axes[2].set_ylabel(r"$\left| S_X^3 \left( f, f' \right) \right|$ / $\left| S_X^3 \left( f_0, f' \right) \right|$",
                       fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

    start_time = time.time()
    for temperature_index, temperature in reverse_enumerate(temperatures):
        print(f"Temperature = {temperature:.4f}")
        current_color = next(colors)
        power_spectrum = get_power_spectrum(algorithm_name, observable_string, output_directory, temperature,
                                            temperature_index, no_of_sites, no_of_sites_string,
                                            external_global_moves_string, no_of_runs, pool, no_of_equilibration_sweeps)
        second_spectrum = get_second_spectrum(algorithm_name, observable_string, output_directory, temperature,
                                              temperature_index, no_of_sites, no_of_sites_string,
                                              external_global_moves_string, no_of_runs, pool,
                                              no_of_equilibration_sweeps)
        power_trispectrum = get_power_trispectrum_nonzero_mode(
            algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
            no_of_sites_string, external_global_moves_string, no_of_runs, pool,
            no_of_trispectrum_auxiliary_frequency_octaves, trispectrum_base_period_shift, target_auxiliary_frequency,
            no_of_equilibration_sweeps)

        # normalise polyspectra with respect to their low-frequency values
        if power_spectrum[1, 0] == 0.0:
            power_spectrum[1] = np.zeros(len(power_spectrum[1]))
        else:
            power_spectrum[1] /= power_spectrum[1, 0]
        if second_spectrum[1, 0] == 0.0:
            second_spectrum[1] = np.zeros(len(second_spectrum[1]))
        else:
            second_spectrum[1] /= second_spectrum[1, 0]
        if power_trispectrum[2][0] == 0.0:
            power_trispectrum[2] = np.zeros(len(power_trispectrum[2]))
        else:
            power_trispectrum[2] /= power_trispectrum[2][0]

        # note frequency indices in order to discard spectrum (trispectrum) for all frequencies >= 0.1 (5.0e-3)  since
        # i) interesting physics at long timescales, and ii) emergent Langevin diffusion breaks down at short timescales
        max_power_spectrum_index = np.argmax(power_spectrum[0] > 1.0e-1) - 1
        max_second_spectrum_index = np.argmax(second_spectrum[0] > 1.0e-1) - 1
        max_power_trispectrum_index = np.argmax(power_trispectrum[1] > 5.0e-3) - 1

        """fit Lorentzian model to power spectrum"""
        final_frequency = 1.0e-2
        (lorentzian_model_parameters, lorentzian_model_errors, lorentzian_model_frequency_values,
         lorentzian_model_spectrum_values) = fit_lorentzian_model_to_spectrum(power_spectrum, final_frequency,
                                                                              no_of_runs)
        axes[0].loglog(power_spectrum[0, :max_power_spectrum_index], power_spectrum[1, :max_power_spectrum_index],
                       color=current_color,
                       label=fr"temperature = {temperature:.4f}; $f_c$ = {lorentzian_model_parameters[1]:.2e} $\pm$ "
                             fr"{lorentzian_model_errors[1]:.2e}; $S_0$ = {lorentzian_model_parameters[0]:.2f} $\pm$ "
                             fr"{lorentzian_model_errors[0]:.2e}")
        axes[0].loglog(lorentzian_model_frequency_values, lorentzian_model_spectrum_values, color='k')

        """fit Lorentzian model to second spectrum"""
        final_frequency = 1.0e-2
        (lorentzian_model_parameters, lorentzian_model_errors, lorentzian_model_frequency_values,
         lorentzian_model_spectrum_values) = fit_lorentzian_model_to_spectrum(second_spectrum, final_frequency,
                                                                              no_of_runs)
        axes[1].loglog(second_spectrum[0, :max_second_spectrum_index], second_spectrum[1, :max_second_spectrum_index],
                       color=current_color,
                       label=fr"temperature = {temperature:.4f}; $f_c$ = {lorentzian_model_parameters[1]:.2e} $\pm$ "
                             fr"{lorentzian_model_errors[1]:.2e}; $S_0$ = {lorentzian_model_parameters[0]:.2f} $\pm$ "
                             fr"{lorentzian_model_errors[0]:.2e}")
        axes[1].loglog(lorentzian_model_frequency_values, lorentzian_model_spectrum_values, color='k')

        """fit 1/f model to trispectrum"""
        (one_over_f_model_parameters, one_over_f_model_errors, one_over_f_model_frequency_values,
         one_over_f_model_spectrum_values) = fit_one_over_f_model_to_trispectrum(power_trispectrum, "upper", 10.0,
                                                                                 no_of_sites, no_of_runs)
        axes[2].loglog(power_trispectrum[1][:max_power_trispectrum_index],
                       power_trispectrum[2][:max_power_trispectrum_index], color=current_color,
                       label=fr"temperature = {temperature:.4f}; f' = {power_trispectrum[0]:.2e}; "
                             fr"$\alpha$ = {one_over_f_model_parameters[1]:.2e} $\pm$ "
                             fr"{one_over_f_model_errors[1]:.2e}")
        axes[2].loglog(one_over_f_model_frequency_values, one_over_f_model_spectrum_values, color='k')

        if temperature_index == 0:
            target_auxiliary_frequency = power_trispectrum[0]
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    figure.tight_layout()
    if observable_string == "cartesian_magnetisation":
        legends = [axes[0].legend(title=r"black lines: fits to $S_X(f)$ / $S_X(f_0) = S_0 / (1 + (f / f_c)^2)$",
                                  loc="lower left", fontsize=8),
                   axes[1].legend(title=r"black lines: fits to $S_X^2(f)$ / $S_X^2(f_0) = S_0 / (1 + (f / f_c)^2)$",
                                  loc="lower left", fontsize=8),
                   axes[2].legend(title=r"black lines: fits to $S_X^3(f, f') \sim f^{-\alpha}$", loc="lower left",
                                  fontsize=8)]
    else:
        legends = [axis.legend(loc="lower left", fontsize=10) for axis in axes]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(1.5) for legend in legends]
    figure.savefig(f"{output_directory}/{observable_string}_polyspectra_{algorithm_name.replace('-', '_')}_"
                   f"{external_global_moves_string}_{no_of_sites_string}_max_trispectrum_shift_eq_"
                   f"{2 ** no_of_trispectrum_auxiliary_frequency_octaves}_x_{trispectrum_base_period_shift}_"
                   f"delta_t.pdf", bbox_inches="tight")
    if no_of_runs > 1:
        pool.close()


def fit_lorentzian_model_to_spectrum(spectrum, final_frequency, no_of_runs):
    initial_frequency = spectrum[0, 0]
    increment = 10.0 ** math.floor(math.log(initial_frequency, 10))
    final_frequency_index = np.argmax(spectrum[0] > final_frequency) - 1
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            # attempt to fit a Lorentzian to the spectrum...
            if no_of_runs > 1:
                parameter_values_and_errors = curve_fit(lorentzian_model, spectrum[0, :final_frequency_index],
                                                        spectrum[1, :final_frequency_index],
                                                        sigma=spectrum[2, :final_frequency_index],
                                                        bounds=(np.array([0.9, initial_frequency]),
                                                                np.array([1.1, final_frequency])))
            else:
                parameter_values_and_errors = curve_fit(lorentzian_model, spectrum[0, :final_frequency_index],
                                                        spectrum[1, :final_frequency_index],
                                                        bounds=(np.array([1.0, initial_frequency]),
                                                                np.array([10.0, final_frequency])))
            parameter_values = parameter_values_and_errors[0]
            parameter_errors = np.sqrt(np.diag(parameter_values_and_errors[1]))
        except (OptimizeWarning, ValueError, RuntimeError) as _:
            # ...but set model parameters to defaults if it fails
            parameter_values = np.array([1.0, 1.0e-2])
            parameter_errors = np.array([0.0, 0.0])
    model_frequency_values = np.arange(initial_frequency, final_frequency, increment)
    model_spectrum_values = lorentzian_model(model_frequency_values, *parameter_values)
    return parameter_values, parameter_errors, model_frequency_values, model_spectrum_values


def lorentzian_model(frequencies, zero_frequency_value, characteristic_frequency):
    return zero_frequency_value / (1.0 + (frequencies / characteristic_frequency) ** 2.0)


def fit_one_over_f_model_to_trispectrum(power_trispectrum, frequency_range, max_model_exponent, no_of_sites,
                                        no_of_runs):
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
            initial_frequency, final_frequency = 2.0e-4, 4.0e-4
        elif no_of_sites <= 16 ** 2:
            initial_frequency, final_frequency = 8.0e-4, 1.5e-3
    else:
        raise SystemExit("Give 'lower' or 'upper' as the value for frequency_range in the "
                         "fit_one_over_f_model_to_trispectrum() method.")
    increment = 10.0 ** math.floor(math.log(initial_frequency, 10)) / 2.0
    initial_frequency_index = np.argmax(power_trispectrum[1] > initial_frequency) - 1
    final_frequency_index = np.argmax(power_trispectrum[1] > final_frequency) - 1
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            # attempt to fit a f^{-alpha} model to the trispectrum...
            if no_of_runs > 32:
                """n.b., no_of_runs > 32 as we found worse results using trispectrum error bars for no_of_runs = 8."""
                parameter_values_and_errors = curve_fit(
                    one_over_f_model, power_trispectrum[1][initial_frequency_index:final_frequency_index],
                    power_trispectrum[2][initial_frequency_index:final_frequency_index],
                    sigma=power_trispectrum[3][initial_frequency_index:final_frequency_index],
                    bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))
            else:
                parameter_values_and_errors = curve_fit(
                    one_over_f_model, power_trispectrum[1][initial_frequency_index:final_frequency_index],
                    power_trispectrum[2][initial_frequency_index:final_frequency_index],
                    bounds=(np.array([0.0, 0.0]), np.array([10.0, max_model_exponent])))
            parameter_values = parameter_values_and_errors[0]
            parameter_errors = np.sqrt(np.diag(parameter_values_and_errors[1]))
        except (OptimizeWarning, ValueError, RuntimeError) as _:
            # ...but set model parameters to defaults if it fails
            parameter_values = np.array([1.0, 0.0])
            parameter_errors = np.array([0.0, 0.0])
    model_frequency_values = np.arange(initial_frequency, final_frequency, increment)
    model_spectrum_values = one_over_f_model(model_frequency_values, *parameter_values)
    return parameter_values, parameter_errors, model_frequency_values, model_spectrum_values


def one_over_f_model(frequencies, scale_factor, exponent):
    return scale_factor / frequencies ** exponent


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 6:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location "
                        "(str) and the observable whose power trispectrum you wish to estimate (str) in the first and "
                        "second positions, respectively.  In addition, you may provide "
                        "no_of_trispectrum_auxiliary_frequency_octaves (default value is 6), "
                        "trispectrum_base_period_shift (default value is 1) and target_auxiliary_frequency (a float "
                        "with default value None) in the third, fourth and fifth positions, respectively.")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file "
              "(str) / the observable whose power trispectrum you wish to estimate (str).  In addition, you may provide"
              " no_of_trispectrum_auxiliary_frequency_octaves (default value is 6), trispectrum_base_period_shift "
              "(default value is 1) and target_auxiliary_frequency (a float with default value None) in the third, "
              "fourth and fifth positions (respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file (str) / the observable whose power trispectrum you wish to estimate (str) / "
              "no_of_trispectrum_auxiliary_frequency_octaves (int).  In addition, you may provide "
              "trispectrum_base_period_shift (default value is 1) and target_auxiliary_frequency (a float with default "
              "value None) in the fourth and fifth positions, respectively.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The first / second / third / fourth must be the location of the "
              "configuration file (str) / the observable whose power trispectrum you wish to estimate (str) / "
              "no_of_trispectrum_auxiliary_frequency_octaves (int) / trispectrum_base_period_shift (int).  In "
              "addition, you may provide target_auxiliary_frequency (a float with default value None) in the fifth "
              "position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    elif len(sys.argv) == 6:
        print("Five positional arguments provided.  The first / second / third / fourth / fifth must be the location "
              "of the configuration file (str) / the observable whose power trispectrum you wish to estimate (str) / "
              "no_of_trispectrum_auxiliary_frequency_octaves (int) / trispectrum_base_period_shift (int) / "
              "target_auxiliary_frequency (float or None).")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]))
