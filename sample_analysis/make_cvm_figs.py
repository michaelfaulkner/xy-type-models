from cramer_von_mises import get_cvm_mag_phase_vs_temperature
from make_twist_figs import compute_physical_time_steps, get_mag_phase_sample_variances_and_errors, \
    get_twist_probabilities_and_errors
from sample_getter import get_no_of_events
from setup_scripts import setup_pool
import importlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(no_of_system_sizes=6):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    linear_system_sizes = [2 ** (index + 2) for index in range(no_of_system_sizes)]
    base_config_file_metrop_low_temps = f"config_files/cvm_figs/4x4_metrop_low_temps.txt"
    base_config_file_metrop_lower_trans = f"config_files/cvm_figs/4x4_metrop_lower_trans.txt"
    base_config_file_metrop_upper_trans = f"config_files/cvm_figs/4x4_metrop_upper_trans.txt"
    base_config_file_metrop_high_temps = f"config_files/cvm_figs/4x4_metrop_high_temps.txt"
    base_config_file_ecmc = f"config_files/cvm_figs/4x4_ecmc.txt"

    """temperatures_ecmc = [*temperatures_metrop_low_temps, *temperatures_metrop_high_temps] in the following lines"""
    (algorithm_name_metrop, sample_directory_4x4_metrop_low_temps, _, _, no_of_equilibration_sweeps_metrop_low_temps,
     no_of_observations_metrop_low_temps, temperatures_metrop_low_temps, _, external_global_moves_string,
     no_of_jobs_metrop_low_temps, _, max_no_of_cpus) = run_script.get_config_data(base_config_file_metrop_low_temps)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_lower_trans, no_of_observations_metrop_lower_trans,
     temperatures_metrop_lower_trans, _, _, no_of_jobs_metrop_lower_trans, _, _) = run_script.get_config_data(
        base_config_file_metrop_lower_trans)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_upper_trans, no_of_observations_metrop_upper_trans,
     temperatures_metrop_upper_trans, _, _, no_of_jobs_metrop_upper_trans, _, _) = run_script.get_config_data(
        base_config_file_metrop_upper_trans)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
     temperatures_metrop_high_temps, _, _, no_of_jobs_metrop_high_temps, _, _) = run_script.get_config_data(
        base_config_file_metrop_high_temps)
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_ecmc, _, _,
     no_of_jobs_ecmc, _, _) = run_script.get_config_data(base_config_file_ecmc)

    output_directory = sample_directory_4x4_metrop_low_temps.replace("/4x4_metrop_low_temps", "")
    sample_directories_metrop_low_temps = [f"{output_directory}/{length}x{length}_metrop_low_temps" for length in
                                           linear_system_sizes]
    sample_directories_metrop_lower_trans = [f"{output_directory}/{length}x{length}_metrop_lower_trans" for length in
                                             linear_system_sizes]
    sample_directories_metrop_upper_trans = [f"{output_directory}/{length}x{length}_metrop_upper_trans" for length in
                                             linear_system_sizes]
    sample_directories_metrop_high_temps = [f"{output_directory}/{length}x{length}_metrop_high_temps" for length in
                                            linear_system_sizes]
    sample_directories_ecmc = [f"{output_directory}/{length}x{length}_ecmc" for length in linear_system_sizes]
    pool = setup_pool(no_of_jobs_metrop_low_temps, max_no_of_cpus)

    temperatures_metrop = [*temperatures_metrop_low_temps, *temperatures_metrop_lower_trans,
                           *temperatures_metrop_upper_trans, *temperatures_metrop_high_temps]
    temperatures_around_transition = temperatures_metrop[3:15]
    approx_transition_temperature = 0.887
    reduced_temperatures_metrop = [temperature / approx_transition_temperature for temperature in temperatures_metrop]
    reduced_temperatures_ecmc = [temperature / approx_transition_temperature for temperature in temperatures_ecmc]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.75), gridspec_kw={'width_ratios': [1.2, 1.0]})
    fig.text(-0.025, 0.9, "d", fontsize=20, weight='bold')
    fig.text(0.525, 0.9, "e", fontsize=20, weight='bold')
    fig.tight_layout(w_pad=3.0)
    axes[0].set_yscale('log')
    axes[0].set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=20, labelpad=-0.5)
    axes[0].set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=20, labelpad=-7.5)
    axes[0].set_ylim([2.0e-4, 2.5e5])
    axes[0].set_yticks([10.0 ** index for index in range(-3, 6)])
    axes[1].set_xlabel(r"$1 / \left( \ln N \right)^2$", fontsize=20, labelpad=4)
    axes[1].set_ylabel(r"$\widetilde{\beta}_{\rm BKT} / \beta_{\rm int}$", fontsize=20, labelpad=-0.5)
    axes[1].set_xlim([0.0, 0.06]), axes[1].set_ylim([0.985, 1.95])
    [axis.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='major', direction='in', width=2.5, length=5, labelsize=18, pad=2.5)
        axis.tick_params(which='minor', direction='in', width=1.5, length=4)

    inset_axis = plt.axes([0.6, 0.55, 0.21, 0.37])
    inset_axis.tick_params(which='major', direction='in', width=2, length=4, labelsize=12, pad=3.5)
    inset_axis.tick_params(which='minor', direction='in', width=0.25, length=3)
    inset_axis.yaxis.set_label_position("right"), inset_axis.yaxis.tick_right()
    inset_axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=14, labelpad=1.0)
    inset_axis.set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=14, labelpad=1.5)
    inset_axis.set_xlim([0.98, 1.54]), inset_axis.set_ylim([0.75, 5.0 * 10 ** 4])
    inset_axis.set_yscale('log')
    inset_axis.set_yticks([1.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4])
    [inset_axis.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"]]

    twinned_axis = axes[1].twiny().twinx()
    twinned_axis.set_xlabel(r"$N$", fontsize=16, labelpad=-0.5, color="red")
    twinned_axis.set_ylabel(r"$n \omega_{\phi_m,n}^2 \left(\beta = \beta_{\rm int} \right)$", fontsize=16,
                            labelpad=-0.5, color="red")
    # twinned_axis.set_yticks([0.0, 0.1, 0.2, 0.3])
    twinned_axis.set_xlim([50.0, 2.25e4]), twinned_axis.set_ylim([1.5, 4.0e4])
    twinned_axis.tick_params(which='major', direction='in', width=2.5, length=5, labelsize=18, pad=2.5, colors="red")
    twinned_axis.tick_params(which='minor', direction='in', width=1.5, length=4, colors="red")
    twinned_axis.tick_params(axis='both', labelcolor='red')
    twinned_axis.spines["right"].set_color("red"), twinned_axis.spines["right"].set_linewidth(3.0)
    twinned_axis.spines["top"].set_color("red"), twinned_axis.spines["top"].set_linewidth(3.0)

    supplementary_fig, supplementary_axis = plt.subplots(1, figsize=(6.25, 4.0))
    supplementary_fig.tight_layout()
    supplementary_axis.set_xlabel(r"$N$", fontsize=20, labelpad=-0.5)
    supplementary_axis.set_ylabel(r"$n \omega_{\phi_m,n}^2 \left(\beta = \beta_{\rm int} \right)$", fontsize=20,
                                  labelpad=-0.5)
    supplementary_axis.tick_params(which='major', direction='in', width=2.5, length=5, labelsize=18, pad=2.5)
    supplementary_axis.tick_params(which='major', axis='x', pad=5)
    supplementary_axis.tick_params(which='minor', direction='in', width=1.5, length=4)
    [supplementary_axis.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"]]

    colors = ["black", "red", "blue", "green", "magenta", "indigo"][:no_of_system_sizes]
    colors.reverse()

    reduced_intersect_temperatures, intersect_values = [], []
    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        print(f"Number of sites, N = {length}x{length}")
        """compute non-used physical time steps for our records"""
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_low_temps[system_size_index], temperatures_metrop_low_temps, length,
            no_of_observations_metrop_low_temps, no_of_jobs_metrop_low_temps, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_lower_trans[system_size_index], temperatures_metrop_lower_trans, length,
            no_of_observations_metrop_lower_trans, no_of_jobs_metrop_lower_trans, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_upper_trans[system_size_index], temperatures_metrop_upper_trans, length,
            no_of_observations_metrop_upper_trans, no_of_jobs_metrop_upper_trans, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_metrop_high_temps, length,
            no_of_observations_metrop_high_temps, no_of_jobs_metrop_high_temps, pool)
        """compute non-used event rates for our records"""
        compute_event_rates(algorithm_name_ecmc, external_global_moves_string, output_directory,
                            sample_directories_ecmc[system_size_index], temperatures_ecmc, length,
                            no_of_observations_ecmc, no_of_jobs_ecmc, pool)

        """compute non-used sample variances for our records"""
        '''_, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_low_temps[system_size_index], temperatures_metrop_low_temps, length,
            no_of_equilibration_sweeps_metrop_low_temps, no_of_observations_metrop_low_temps,
            no_of_jobs_metrop_low_temps, pool))
        _, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_lower_trans[system_size_index], temperatures_metrop_lower_trans, length,
            no_of_equilibration_sweeps_metrop_lower_trans, no_of_observations_metrop_lower_trans,
            no_of_jobs_metrop_lower_trans, pool))
        _, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_upper_trans[system_size_index], temperatures_metrop_upper_trans, length,
            no_of_equilibration_sweeps_metrop_upper_trans, no_of_observations_metrop_upper_trans,
            no_of_jobs_metrop_upper_trans, pool))
        _, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_metrop_high_temps, length,
            no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
            no_of_jobs_metrop_high_temps, pool))
        _, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name_ecmc, external_global_moves_string, output_directory,
            sample_directories_ecmc[system_size_index], temperatures_ecmc, length, no_of_equilibration_sweeps_ecmc,
            no_of_observations_ecmc, no_of_jobs_ecmc, pool))'''

        """compute non-used (in this script) twist probabilities for later use in make_twist_figs -- this optimises 
            our usage of the scratch space on BlueCrystal 4, ACRC, University of Bristol"""
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
            temperatures_metrop_low_temps, length, no_of_observations_metrop_low_temps, no_of_jobs_metrop_low_temps,
            pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_lower_trans[system_size_index],
            temperatures_metrop_lower_trans, length, no_of_observations_metrop_lower_trans,
            no_of_jobs_metrop_lower_trans, pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_upper_trans[system_size_index],
            temperatures_metrop_upper_trans, length, no_of_observations_metrop_upper_trans,
            no_of_jobs_metrop_upper_trans, pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            temperatures_metrop_high_temps, length, no_of_observations_metrop_high_temps, no_of_jobs_metrop_high_temps,
            pool)

        cvms_metrop_low_temps, cvm_errors_metrop_low_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_low_temps, no_of_observations_metrop_low_temps,
            temperatures_metrop_low_temps, external_global_moves_string, no_of_jobs_metrop_low_temps, pool, length)
        cvms_metrop_lower_trans, cvm_errors_metrop_lower_trans = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_lower_trans[system_size_index],
            no_of_equilibration_sweeps_metrop_lower_trans, no_of_observations_metrop_lower_trans,
            temperatures_metrop_lower_trans, external_global_moves_string, no_of_jobs_metrop_lower_trans, pool, length)
        cvms_metrop_upper_trans, cvm_errors_metrop_upper_trans = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_upper_trans[system_size_index],
            no_of_equilibration_sweeps_metrop_upper_trans, no_of_observations_metrop_upper_trans,
            temperatures_metrop_upper_trans, external_global_moves_string, no_of_jobs_metrop_upper_trans, pool, length)
        cvms_metrop_high_temps, cvm_errors_metrop_high_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
            temperatures_metrop_high_temps, external_global_moves_string, no_of_jobs_metrop_high_temps, pool, length)
        cvms_ecmc, cvm_errors_ecmc = get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, sample_directories_ecmc[system_size_index],
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_ecmc,
            external_global_moves_string, no_of_jobs_ecmc, pool, length)

        cvms_metrop = np.array([*cvms_metrop_low_temps, *cvms_metrop_lower_trans, *cvms_metrop_upper_trans,
                                *cvms_metrop_high_temps])
        cvm_errors_metrop = np.array([*cvm_errors_metrop_low_temps, *cvm_errors_metrop_lower_trans,
                                      *cvm_errors_metrop_upper_trans, *cvm_errors_metrop_high_temps])

        """estimate intercept temperatures - we fit a 2nd-order polynomial to log(1 / CvM) as this was the smoothest 
            fitting we found - ***NOTE THAT*** using three temperature values with a 2nd-order polynomial implies no 
            errors on the fitting parameters"""
        cvms_around_transition = cvms_metrop[3:15]
        cvm_errors_around_transition = cvm_errors_metrop[3:15]
        if system_size_index > 0:
            current_fitting_cvms = cvms_around_transition[lower_fitting_index:upper_fitting_index]
            current_fitting_cvm_errors = cvm_errors_around_transition[lower_fitting_index:upper_fitting_index]
            current_polynomial_fit = np.poly1d(np.polyfit(fitting_temps, np.log(current_fitting_cvms), 2,
                                                          w=current_fitting_cvm_errors / np.abs(current_fitting_cvms)))
            continuous_temperatures = np.linspace(fitting_temps[0], fitting_temps[-1], 100)
            """the following commented-out code plots the fittings on top of the data in the inset figure"""
            '''inset_axis.plot(continuous_temperatures / approx_transition_temperature,
                            np.exp(previous_polynomial_fit(continuous_temperatures)), 
                            color=colors[system_size_index - 1], linestyle="-")
            inset_axis.plot(continuous_temperatures / approx_transition_temperature,
                            np.exp(current_polynomial_fit(continuous_temperatures)), color=colors[system_size_index],
                            linestyle="-")'''
            intersect_index = np.argwhere(np.diff(
                np.sign(previous_polynomial_fit(continuous_temperatures) -
                        current_polynomial_fit(continuous_temperatures)))).flatten()[0]
            intersect_temperature = continuous_temperatures[intersect_index]
            intersect_value = np.exp(current_polynomial_fit(continuous_temperatures[intersect_index]))
            inset_axis.vlines(intersect_temperature / approx_transition_temperature, intersect_value, 5.0e5,
                              colors="gray", linestyles='solid')
            reduced_intersect_temperatures.append(intersect_temperature / approx_transition_temperature)
            intersect_values.append(intersect_value)

        if system_size_index == 0:
            lower_fitting_index, upper_fitting_index = 9, 12
        if system_size_index == 1:
            lower_fitting_index, upper_fitting_index = 5, 8
        if system_size_index == 2:
            lower_fitting_index, upper_fitting_index = 2, 5
        if system_size_index == 3:
            lower_fitting_index, upper_fitting_index = 1, 4
        if system_size_index == 4:
            lower_fitting_index, upper_fitting_index = 0, 3
        fitting_temps = temperatures_around_transition[lower_fitting_index:upper_fitting_index]
        previous_fitting_cvms = cvms_around_transition[lower_fitting_index:upper_fitting_index]
        previous_fitting_cvm_errors = cvm_errors_around_transition[lower_fitting_index:upper_fitting_index]
        previous_polynomial_fit = np.poly1d(np.polyfit(fitting_temps, np.log(previous_fitting_cvms), 2,
                                                       w=previous_fitting_cvm_errors / np.abs(previous_fitting_cvms)))
        """the following commented-out code was an attempt at using scipy.optimize.curve_fit()"""
        '''popt, pcov = curve_fit(quadratic_fit, fitting_temps, np.log(previous_fitting_cvms),
                               sigma=previous_fitting_cvm_errors / previous_fitting_cvms)'''

        """plot CvM data"""
        if system_size_index == 5:
            """lowest-temperature results did not converge for L = 128"""
            axes[0].errorbar(reduced_temperatures_metrop[1:], cvms_metrop[1:], cvm_errors_metrop[1:], marker=".",
                             markersize=10, color=colors[system_size_index], linestyle="None",
                             label=fr"$N$ = {length}x{length} Metrop.")
        else:
            axes[0].errorbar(reduced_temperatures_metrop, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                             color=colors[system_size_index], linestyle="None",
                             label=fr"$N$ = {length}x{length} Metrop.")
        inset_axis.errorbar(reduced_temperatures_metrop, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                            color=colors[system_size_index], linestyle="None")
        axes[0].errorbar(reduced_temperatures_ecmc, cvms_ecmc, cvm_errors_ecmc, marker="*", markersize=8,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length} ECMC")

    """fit reduced intercept temperatures"""
    inverse_log_squared_system_sizes = [1.0 / np.log(length ** 2) ** 2 for index, length in
                                        enumerate(linear_system_sizes) if index > 0]
    inverse_log_squared_system_sizes.reverse(), reduced_intersect_temperatures.reverse()
    continuous_inverse_log_squared_system_sizes = np.linspace(0.0, inverse_log_squared_system_sizes[-1] + 0.1, 100)
    temp_polyfit_params, temp_polyfit_cov = np.polyfit(inverse_log_squared_system_sizes, reduced_intersect_temperatures,
                                                       1, cov=True)
    temp_polyfit_errors = np.sqrt(np.diag(temp_polyfit_cov))
    print(f"Thermodynamic reduced intercept temperature = {temp_polyfit_params[1]} +- {temp_polyfit_errors[1]}")
    left_legend_1 = axes[1].plot(inverse_log_squared_system_sizes, reduced_intersect_temperatures, marker=".",
                                 markersize=10, color="black", linestyle="None", label="intercept temps")
    temp_polynomial_fit = np.poly1d(temp_polyfit_params)
    left_legend_2 = axes[1].plot(continuous_inverse_log_squared_system_sizes,
                                 temp_polynomial_fit(continuous_inverse_log_squared_system_sizes), linestyle="--",
                                 color="black", label="temp. fit")

    """fit intercept values to power law"""
    system_sizes = [length ** 2 for index, length in enumerate(linear_system_sizes) if index > 0]
    continuous_system_sizes = np.linspace(0.5 * system_sizes[0], 2.0 * system_sizes[-1], 100)
    value_polyfit_params, value_polyfit_cov = np.polyfit(np.log(system_sizes), np.log(intersect_values), 1, cov=True)
    value_polyfit_errors = np.sqrt(np.diag(value_polyfit_cov))
    print(f"Intercept-value power law = {np.exp(value_polyfit_params[1])} x N^a, where a = {value_polyfit_params[0]} "
          f"+- {value_polyfit_errors[0]}")
    '''right_legend_1 = twinned_axis.loglog(system_sizes, intersect_values, marker="*", markersize=10, color="red",
                                         linestyle="None", label="intercept values")
    right_legend_2 = axes[1].loglog(continuous_system_sizes, np.exp(value_polyfit_params[1]) *
                                    continuous_system_sizes ** value_polyfit_params[0], linestyle="--", color="red",
                                    label=rf"\sim N^{value_polyfit_params[0]:.3f}")'''

    supplementary_legend_1 = supplementary_axis.loglog(system_sizes, intersect_values, marker=".", markersize=10,
                                                       color="black", linestyle="None", label="intercept values")
    supplementary_legend_2 = supplementary_axis.loglog(
        continuous_system_sizes, np.exp(value_polyfit_params[1]) * continuous_system_sizes ** value_polyfit_params[0],
        linestyle="--", color="black", label=fr"$\sim N^\alpha$, $\alpha = {value_polyfit_params[0]:.3f}$")

    right_legend_1 = twinned_axis.loglog(system_sizes, intersect_values, marker="*", markersize=10, color="red",
                                         linestyle="None", label="intercept values")
    right_legend_2 = twinned_axis.loglog(
        continuous_system_sizes, np.exp(value_polyfit_params[1]) * continuous_system_sizes ** value_polyfit_params[0],
        linestyle="--", color="red", label=fr"$\sim N^\alpha$, $\alpha = {value_polyfit_params[0]:.3f}$")

    """output reduced intercept temperatures, inverse intercept values & thermodynamic reduced intercept temperature"""
    intercept_temps_file = open(f"{output_directory}/cramer_von_mises_intercept_temps_and_inverse_values_xy_gaussian_"
                                f"noise_metropolis_w_twists.tsv", "w")
    intercept_temps_file.write("# system size, N".ljust(30) + "reduced intercept temp.".ljust(40) +
                               "inverse intercept value" + "\n")
    linear_system_sizes.reverse()
    linear_system_sizes.pop()
    linear_system_sizes.reverse(), reduced_intersect_temperatures.reverse(), intersect_values.reverse()
    for system_size_index, length in enumerate(linear_system_sizes):
        intercept_temps_file.write(f"{length}x{length}".ljust(30) +
                                   f"{reduced_intersect_temperatures[system_size_index]:.14e}".ljust(40) +
                                   f"{intersect_values[system_size_index]:.14e}" + "\n")
    intercept_temps_file.write("\n\n" + f"# Thermodynamic reduced intercept temperature = {temp_polyfit_params[1]} +- "
                                        f"{temp_polyfit_errors[1]}")
    intercept_temps_file.write("\n" + f"# Intercept-value power law = {np.exp(value_polyfit_params[1])} x N^a, "
                                      f"where a = {value_polyfit_params[0]} +- {value_polyfit_errors[0]}")
    intercept_temps_file.close()

    combined_legends_data = left_legend_1 + left_legend_2 + right_legend_1 + right_legend_2
    combined_legends_labels = [legend.get_label() for legend in combined_legends_data]
    legends = [axes[0].legend(loc="upper right", ncol=2, fontsize=7.625),
               axes[1].legend(combined_legends_data, combined_legends_labels, loc=[0.625, 0.05], fontsize=8.5)]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(3.0) for legend in legends]
    fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc.pdf",
                bbox_inches="tight")

    combined_supplementary_legends_data = supplementary_legend_1 + supplementary_legend_2
    combined_supplementary_legends_labels = [legend.get_label() for legend in combined_supplementary_legends_data]
    supplementary_legend = supplementary_axis.legend(
        combined_supplementary_legends_data, combined_supplementary_legends_labels, loc="lower right", fontsize=12)
    supplementary_legend.get_frame().set_edgecolor("k"), supplementary_legend.get_frame().set_lw(3.0)
    supplementary_fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_"
                              f"intercept_values.pdf", bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


def compute_event_rates(algorithm_name, external_global_moves_string, output_directory, sample_directory, temperatures,
                        length, no_of_observations, no_of_jobs, pool):
    try:
        with open(f"{output_directory}/event_rates_{algorithm_name.replace('-', '_')}_{external_global_moves_string}_"
                  f"{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_"
                  f"{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "r") as _:
            pass
    except IOError:
        event_rates_file = open(
            f"{output_directory}/event_rates_{algorithm_name.replace('-', '_')}_{external_global_moves_string}_"
            f"{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_{no_of_observations}_"
            f"obs_{no_of_jobs}_jobs.tsv", "w")
        event_rates_file.write("# temperature".ljust(30) + "event rate".ljust(30) + "event-rate error" + "\n")
        for temperature_index, temperature in enumerate(temperatures):
            event_rate_vs_job = np.array(pool.starmap(get_no_of_events,
                                                      [(f"{sample_directory}/job_{job_index}", temperature_index)
                                                       for job_index in range(no_of_jobs)])) / no_of_observations
            event_rates_file.write(f"{temperature:.14e}".ljust(30) + f"{np.mean(event_rate_vs_job):.14e}".ljust(30) +
                                   f"{np.std(event_rate_vs_job) / len(event_rate_vs_job) ** 0.5:.14e}" + "\n")
        event_rates_file.close()


def power_law_fit(one_dimensional_array, alpha, beta):
    return alpha * one_dimensional_array ** beta


def quadratic_fit(one_dimensional_array, a, b, c):
    return a * one_dimensional_array ** 2 + b * one_dimensional_array + c


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional arguments permitted.  None are required but you may "
                        "provide no_of_system_sizes, which must be an integer greater than 0 and less than 8 (default "
                        "value is 6).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be no_of_system_sizes - which must be an integer greater "
              "than 0 and less than 8 (default value is 6).")
        chosen_no_of_system_sizes = int(sys.argv[1])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 7:
            raise Exception("InterfaceError: no_of_system_sizes must be an integer greater than 0 and less than 8 "
                            "(default value is 6).")
        main(chosen_no_of_system_sizes)
    else:
        print("No positional arguments provided.  None are required but you may provide no_of_system_sizes, which must "
              "be an integer greater than 0 and less than 8 (default value is 6).")
        main()
