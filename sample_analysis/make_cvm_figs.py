from cramer_von_mises import get_cvm_mag_phase_vs_temperature
from make_twist_figs import compute_physical_time_steps, get_mag_phase_simulation_variances_and_errors, \
    get_twist_probabilities_and_errors
from sample_getter import get_event_rate
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


def main(symmetry_breaking_paper=True, no_of_system_sizes=6):
    """For 'Symmetry breaking at a topological phase transition', set symmetry_breaking_paper=True;
        for 'Emergent electrostatics in...' paper, set symmetry_breaking_paper=False."""
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    linear_system_sizes = [2 ** (index + 2) for index in range(no_of_system_sizes)]
    base_config_file_metrop_low_temps = f"config_files/cvm_figs/4x4_metrop_low_temps.txt"
    base_config_file_metrop_lower_trans = f"config_files/cvm_figs/4x4_metrop_lower_trans.txt"
    base_config_file_metrop_upper_trans = f"config_files/cvm_figs/4x4_metrop_upper_trans.txt"
    base_config_file_metrop_high_temps = f"config_files/cvm_figs/4x4_metrop_high_temps.txt"
    base_config_file_ecmc = f"config_files/cvm_figs/4x4_ecmc.txt"

    """temperatures_ecmc = [*temperatures_metrop_low_temps, *temperatures_metrop_high_temps] in the following lines"""
    (algorithm_name_metrop, sample_directory_4x4_metrop_low_temps, _, _, no_of_equilibration_sweeps_metrop_low_temps,
     no_of_samples_metrop_small_systems_low_temps, temperatures_metrop_low_temps, _, external_global_moves_string,
     no_of_runs_metrop_low_temps, _, max_no_of_cpus) = run_script.get_config_data(base_config_file_metrop_low_temps)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_lower_trans, no_of_samples_metrop_lower_trans,
     temperatures_metrop_lower_trans, _, _, no_of_runs_metrop_lower_trans, _, _) = run_script.get_config_data(
        base_config_file_metrop_lower_trans)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_upper_trans, no_of_samples_metrop_upper_trans,
     temperatures_metrop_upper_trans, _, _, no_of_runs_metrop_upper_trans, _, _) = run_script.get_config_data(
        base_config_file_metrop_upper_trans)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_high_temps, no_of_samples_metrop_high_temps,
     temperatures_metrop_high_temps, _, _, no_of_runs_metrop_high_temps, _, _) = run_script.get_config_data(
        base_config_file_metrop_high_temps)
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_samples_ecmc, temperatures_ecmc, _, _,
     no_of_runs_ecmc, _, _) = run_script.get_config_data(base_config_file_ecmc)

    """We also define an additional 'lowest-temps' case for L = 128.  This is because we split the lowest four 
        temperatures into two directories in this case.  This keeps within the two-week time limit on BlueCrystal 4, as 
        this larger system size uses more CPU time at fixed simulation timescale, but also requires longer simulation 
        timescales for CvM convergence.  Due to these longer required simulation timescales, L = 64 and 128 both 
        require the additional field no_of_samples_metrop_large_systems_low_temps."""
    base_config_file_128x128_metrop_low_temps_lower = f"config_files/cvm_figs/128x128_metrop_lowest_temps.txt"
    base_config_file_128x128_metrop_low_temps_upper = f"config_files/cvm_figs/128x128_metrop_low_temps.txt"
    temperatures_128x128_metrop_low_temps_lower = run_script.get_config_data(
        base_config_file_128x128_metrop_low_temps_lower)[6]
    (_, _, _, _, _, no_of_samples_metrop_large_systems_low_temps, temperatures_128x128_metrop_low_temps_upper, _,
     _, _, _, _) = run_script.get_config_data(base_config_file_128x128_metrop_low_temps_upper)
    no_of_samples_metrop_low_temps = [no_of_samples_metrop_small_systems_low_temps if index < 4 else
                                           no_of_samples_metrop_large_systems_low_temps for index in
                                           range(len(linear_system_sizes))]

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
    """The following line defines an additional 'lowest-temps' sample directory for L = 128.  This is because we split 
        the lowest four temperatures into two directories in this case.  This keeps within the two-week time limit on 
        BlueCrystal 4, as this larger system size uses more CPU time at fixed simulation timescale, but also requires 
        longer simulation timescales for CvM convergence."""
    sample_directory_128x128_metrop_low_temps_lower = f"{output_directory}/128x128_metrop_lowest_temps"

    pool = setup_pool(no_of_runs_metrop_low_temps, max_no_of_cpus)
    temperatures_metrop = [*temperatures_metrop_low_temps, *temperatures_metrop_lower_trans,
                           *temperatures_metrop_upper_trans, *temperatures_metrop_high_temps]
    temperatures_around_transition = temperatures_metrop[3:15]
    approx_transition_temperature = 0.887
    reduced_temperatures_metrop = [temperature / approx_transition_temperature for temperature in temperatures_metrop]
    reduced_temperatures_ecmc = [temperature / approx_transition_temperature for temperature in temperatures_ecmc]

    if symmetry_breaking_paper:
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.85), gridspec_kw={'width_ratios': [1.2, 1.0]})
        fig.text(-0.045, 0.8975, "(a)", fontsize=20)
        fig.text(0.51, 0.8975, "(b)", fontsize=20)
        fig.tight_layout(w_pad=3.0)
        main_axis = axes[0]
        main_axis_font_size = 20.0
        main_axis_tick_label_size = 18.0
        main_axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=main_axis_font_size, labelpad=-0.5)
        axes[1].set_xlabel(r"$1 / \left( \ln N \right)^2$", fontsize=20, labelpad=4)
        axes[1].set_ylabel(r"$\widetilde{\beta}_{\rm BKT} / \beta_{\rm int}$", fontsize=20, labelpad=-0.5)
        axes[1].set_xlim([0.0, 0.06]), axes[1].set_ylim([0.925, 1.65])
    else:
        fig, axis = plt.subplots(1, figsize=(5.45454545, 4.5))
        fig.tight_layout()
        main_axis = axis
        axes = [axis]
        main_axis_font_size = 25.0
        main_axis_tick_label_size = 20.0
        main_axis.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        main_axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta$", fontsize=main_axis_font_size,
                             labelpad=-0.5)
    main_axis.set_yscale('log')
    main_axis.set_ylabel(r"$\langle n \omega_{\phi_{\mathbf{m}},n}^2 \rangle$", fontsize=main_axis_font_size,
                         labelpad=-7.5)
    main_axis.set_ylim([2.0e-4, 2.5e5])
    main_axis.set_yticks([10.0 ** index for index in range(-3, 6)])
    [axis.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='major', direction='in', width=2.5, length=5, labelsize=main_axis_tick_label_size,
                         pad=2.5)
        axis.tick_params(which='minor', direction='in', width=1.5, length=4)

    if symmetry_breaking_paper:
        inset_axis_1 = plt.axes((0.28725, 0.5675, 0.2125, 0.3625))
        inset_axis_1_font_size = 12.5
        inset_axis_1_tick_label_size = 12.0
        inset_axis_1.set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=inset_axis_1_font_size, labelpad=1.0)
    else:
        inset_axis_1 = plt.axes((0.52, 0.55, 0.42, 0.39))
        inset_axis_1_font_size = 15.0
        inset_axis_1_tick_label_size = 15.0
        inset_axis_1.set_xlabel(r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta$", fontsize=inset_axis_1_font_size,
                                labelpad=1.0)
    inset_axis_1.tick_params(which='major', direction='in', width=2, length=4, labelsize=inset_axis_1_tick_label_size,
                             pad=3.5)
    inset_axis_1.tick_params(which='minor', direction='in', width=0.25, length=3)
    inset_axis_1.set_ylabel(r"$\langle n \omega_{\phi_{\mathbf{m}},n}^2 \rangle$", fontsize=inset_axis_1_font_size,
                            labelpad=1.5)
    inset_axis_1.set_xlim([0.975, 1.535]), inset_axis_1.set_ylim([0.9, 5.0 * 10 ** 4])
    inset_axis_1.set_xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
    inset_axis_1.set_yscale('log')
    inset_axis_1.set_yticks([1.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4])
    [inset_axis_1.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"]]

    if symmetry_breaking_paper:
        inset_axis_2 = plt.axes((0.6, 0.58, 0.2, 0.35))
        inset_axis_2.tick_params(which='major', direction='in', width=2, length=4, labelsize=12, pad=3.5)
        inset_axis_2.tick_params(which='minor', direction='in', width=0.25, length=3)
        inset_axis_2.yaxis.set_label_position("right"), inset_axis_2.yaxis.tick_right()
        inset_axis_2.set_xlabel(r"$N$", fontsize=14, labelpad=1.0)
        inset_axis_2.set_ylabel(r"$\langle n \omega_{\phi_\mathbf{m},n}^2 \left(\beta = \beta_{\rm int} \right) \rangle$",
                                fontsize=11.5, labelpad=1.5)
        [inset_axis_2.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"]]
    else:
        inset_axis_2 = None

    supplementary_fig, supplementary_axis = plt.subplots(1, figsize=(6.25, 4.0))
    supplementary_fig.tight_layout()
    supplementary_axis.set_xlabel(r"$N$", fontsize=20, labelpad=-0.5)
    supplementary_axis.set_ylabel(r"$\langle n \omega_{\phi_\mathbf{m},n}^2 \left(\beta = \beta_{\rm int} \right) "
                                  r"\rangle$", fontsize=20, labelpad=-0.5)
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
        if length < 128:
            compute_physical_time_steps(
                algorithm_name_metrop, external_global_moves_string, output_directory,
                sample_directories_metrop_low_temps[system_size_index], temperatures_metrop_low_temps, length,
                no_of_samples_metrop_low_temps[system_size_index], no_of_runs_metrop_low_temps, pool)
        else:
            compute_physical_time_steps(
                algorithm_name_metrop, external_global_moves_string, output_directory,
                sample_directory_128x128_metrop_low_temps_lower, temperatures_128x128_metrop_low_temps_lower,
                length, no_of_samples_metrop_low_temps[system_size_index], no_of_runs_metrop_low_temps, pool)
            compute_physical_time_steps(
                algorithm_name_metrop, external_global_moves_string, output_directory,
                sample_directories_metrop_low_temps[system_size_index], temperatures_128x128_metrop_low_temps_upper,
                length, no_of_samples_metrop_low_temps[system_size_index], no_of_runs_metrop_low_temps, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_lower_trans[system_size_index], temperatures_metrop_lower_trans, length,
            no_of_samples_metrop_lower_trans, no_of_runs_metrop_lower_trans, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_upper_trans[system_size_index], temperatures_metrop_upper_trans, length,
            no_of_samples_metrop_upper_trans, no_of_runs_metrop_upper_trans, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_metrop_high_temps, length,
            no_of_samples_metrop_high_temps, no_of_runs_metrop_high_temps, pool)

        """compute non-used event rates for our records"""
        compute_event_rates(algorithm_name_ecmc, external_global_moves_string, output_directory,
                            sample_directories_ecmc[system_size_index], temperatures_ecmc, length,
                            no_of_samples_ecmc, no_of_runs_ecmc, pool)

        """compute non-used simulation variances for our records"""
        if length < 128:
            _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
                algorithm_name_metrop, external_global_moves_string, output_directory,
                sample_directories_metrop_low_temps[system_size_index], temperatures_metrop_low_temps, length,
                no_of_equilibration_sweeps_metrop_low_temps, no_of_samples_metrop_low_temps[system_size_index],
                no_of_runs_metrop_low_temps, pool))
        else:
            _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
                algorithm_name_metrop, external_global_moves_string, output_directory,
                sample_directory_128x128_metrop_low_temps_lower, temperatures_128x128_metrop_low_temps_lower, length,
                no_of_equilibration_sweeps_metrop_low_temps, no_of_samples_metrop_low_temps[system_size_index],
                no_of_runs_metrop_low_temps, pool))
            _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
                algorithm_name_metrop, external_global_moves_string, output_directory,
                sample_directories_metrop_low_temps[system_size_index], temperatures_128x128_metrop_low_temps_upper,
                length, no_of_equilibration_sweeps_metrop_low_temps,
                no_of_samples_metrop_low_temps[system_size_index], no_of_runs_metrop_low_temps, pool))
        _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_lower_trans[system_size_index], temperatures_metrop_lower_trans, length,
            no_of_equilibration_sweeps_metrop_lower_trans, no_of_samples_metrop_lower_trans,
            no_of_runs_metrop_lower_trans, pool))
        _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_upper_trans[system_size_index], temperatures_metrop_upper_trans, length,
            no_of_equilibration_sweeps_metrop_upper_trans, no_of_samples_metrop_upper_trans,
            no_of_runs_metrop_upper_trans, pool))
        _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_metrop_high_temps, length,
            no_of_equilibration_sweeps_metrop_high_temps, no_of_samples_metrop_high_temps,
            no_of_runs_metrop_high_temps, pool))
        _, _ = np.array(get_mag_phase_simulation_variances_and_errors(
            algorithm_name_ecmc, external_global_moves_string, output_directory,
            sample_directories_ecmc[system_size_index], temperatures_ecmc, length, no_of_equilibration_sweeps_ecmc,
            no_of_samples_ecmc, no_of_runs_ecmc, pool))

        """compute non-used (in this script) twist probabilities for later use in make_twist_figs -- this optimises 
            our usage of the scratch space on BlueCrystal 4, ACRC, University of Bristol"""
        if length < 128:
            _, _ = get_twist_probabilities_and_errors(
                algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
                temperatures_metrop_low_temps, length, no_of_samples_metrop_low_temps[system_size_index],
                no_of_runs_metrop_low_temps, pool)
        else:
            _, _ = get_twist_probabilities_and_errors(
                algorithm_name_metrop, output_directory, sample_directory_128x128_metrop_low_temps_lower,
                temperatures_128x128_metrop_low_temps_lower, length,
                no_of_samples_metrop_low_temps[system_size_index], no_of_runs_metrop_low_temps, pool)
            _, _ = get_twist_probabilities_and_errors(
                algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
                temperatures_128x128_metrop_low_temps_upper, length,
                no_of_samples_metrop_low_temps[system_size_index], no_of_runs_metrop_low_temps, pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_lower_trans[system_size_index],
            temperatures_metrop_lower_trans, length, no_of_samples_metrop_lower_trans,
            no_of_runs_metrop_lower_trans, pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_upper_trans[system_size_index],
            temperatures_metrop_upper_trans, length, no_of_samples_metrop_upper_trans,
            no_of_runs_metrop_upper_trans, pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            temperatures_metrop_high_temps, length, no_of_samples_metrop_high_temps, no_of_runs_metrop_high_temps,
            pool)

        if length < 128:
            cvms_metrop_low_temps, cvm_errors_metrop_low_temps = get_cvm_mag_phase_vs_temperature(
                algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
                no_of_equilibration_sweeps_metrop_low_temps, no_of_samples_metrop_low_temps[system_size_index],
                temperatures_metrop_low_temps, external_global_moves_string, no_of_runs_metrop_low_temps, pool, length)
        else:
            cvms_128x128_metrop_lowest_temps, cvm_errors_128x128_metrop_lowest_temps = get_cvm_mag_phase_vs_temperature(
                algorithm_name_metrop, output_directory, sample_directory_128x128_metrop_low_temps_lower,
                no_of_equilibration_sweeps_metrop_low_temps, no_of_samples_metrop_low_temps[system_size_index],
                temperatures_128x128_metrop_low_temps_lower, external_global_moves_string, no_of_runs_metrop_low_temps,
                pool, length)
            cvms_128x128_metrop_low_temps, cvm_errors_128x128_metrop_low_temps = get_cvm_mag_phase_vs_temperature(
                algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
                no_of_equilibration_sweeps_metrop_low_temps, no_of_samples_metrop_low_temps[system_size_index],
                temperatures_128x128_metrop_low_temps_upper, external_global_moves_string, no_of_runs_metrop_low_temps,
                pool, length)
            cvms_metrop_low_temps = [*cvms_128x128_metrop_lowest_temps, *cvms_128x128_metrop_low_temps]
            cvm_errors_metrop_low_temps = [*cvm_errors_128x128_metrop_lowest_temps,
                                           *cvm_errors_128x128_metrop_low_temps]
        cvms_metrop_lower_trans, cvm_errors_metrop_lower_trans = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_lower_trans[system_size_index],
            no_of_equilibration_sweeps_metrop_lower_trans, no_of_samples_metrop_lower_trans,
            temperatures_metrop_lower_trans, external_global_moves_string, no_of_runs_metrop_lower_trans, pool, length)
        cvms_metrop_upper_trans, cvm_errors_metrop_upper_trans = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_upper_trans[system_size_index],
            no_of_equilibration_sweeps_metrop_upper_trans, no_of_samples_metrop_upper_trans,
            temperatures_metrop_upper_trans, external_global_moves_string, no_of_runs_metrop_upper_trans, pool, length)
        cvms_metrop_high_temps, cvm_errors_metrop_high_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_high_temps, no_of_samples_metrop_high_temps,
            temperatures_metrop_high_temps, external_global_moves_string, no_of_runs_metrop_high_temps, pool, length)
        cvms_ecmc, cvm_errors_ecmc = get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, sample_directories_ecmc[system_size_index],
            no_of_equilibration_sweeps_ecmc, no_of_samples_ecmc, temperatures_ecmc,
            external_global_moves_string, no_of_runs_ecmc, pool, length)

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
            '''inset_axis_1.plot(continuous_temperatures / approx_transition_temperature,
                            np.exp(previous_polynomial_fit(continuous_temperatures)),
                            color=colors[system_size_index - 1], linestyle="-")
            inset_axis_1.plot(continuous_temperatures / approx_transition_temperature,
                            np.exp(current_polynomial_fit(continuous_temperatures)), color=colors[system_size_index],
                            linestyle="-")'''
            intersect_index = np.argwhere(np.diff(
                np.sign(previous_polynomial_fit(continuous_temperatures) -
                        current_polynomial_fit(continuous_temperatures)))).flatten()[0]
            intersect_temperature = continuous_temperatures[intersect_index]
            intersect_value = np.exp(current_polynomial_fit(continuous_temperatures[intersect_index]))
            inset_axis_1.vlines(intersect_temperature / approx_transition_temperature, intersect_value, 5.0e5,
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
            lower_fitting_index, upper_fitting_index = 1, 4
        fitting_temps = temperatures_around_transition[lower_fitting_index:upper_fitting_index]
        previous_fitting_cvms = cvms_around_transition[lower_fitting_index:upper_fitting_index]
        previous_fitting_cvm_errors = cvm_errors_around_transition[lower_fitting_index:upper_fitting_index]
        previous_polynomial_fit = np.poly1d(np.polyfit(fitting_temps, np.log(previous_fitting_cvms), 2,
                                                       w=previous_fitting_cvm_errors / np.abs(previous_fitting_cvms)))

        """plot CvM data"""
        if symmetry_breaking_paper:
            metrop_legend_label = fr"$N$ = {length}x{length} Metrop."
        else:
            metrop_legend_label = fr"$N$ = {length}x{length}"

        if system_size_index == 5:
            """lowest-temperature results did not converge for L = 128"""
            main_axis.errorbar(reduced_temperatures_metrop[1:], cvms_metrop[1:], cvm_errors_metrop[1:], marker=".",
                               markersize=10, color=colors[system_size_index], linestyle="None",
                               label=metrop_legend_label)
        else:
            main_axis.errorbar(reduced_temperatures_metrop, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                               color=colors[system_size_index], linestyle="None", label=metrop_legend_label)
        inset_axis_1.errorbar(reduced_temperatures_metrop, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                              color=colors[system_size_index], linestyle="None")
        if symmetry_breaking_paper:
            main_axis.errorbar(reduced_temperatures_ecmc, cvms_ecmc, cvm_errors_ecmc, marker="*", markersize=8,
                               color=colors[system_size_index], linestyle="None",
                               label=fr"$N$ = {length}x{length} ECMC")
        else:
            main_axis.errorbar(reduced_temperatures_ecmc, cvms_ecmc, cvm_errors_ecmc, marker="*", markersize=8,
                               color=colors[system_size_index], linestyle="None")

    if symmetry_breaking_paper:
        """fit reduced intercept temperatures..."""
        inverse_squared_log_system_sizes = [1.0 / np.log(length ** 2) ** 2 for index, length in
                                            enumerate(linear_system_sizes) if index > 0]
        inverse_squared_log_system_sizes.reverse(), reduced_intersect_temperatures.reverse()
        continuous_inverse_squared_log_system_sizes = np.linspace(0.0, inverse_squared_log_system_sizes[-1] + 0.1, 100)
        temp_polyfit_params, temp_polyfit_cov = np.polyfit(
            inverse_squared_log_system_sizes, reduced_intersect_temperatures, 1, cov=True)
        temp_polyfit_errors = np.sqrt(np.diag(temp_polyfit_cov))
        temp_polynomial_fit = np.poly1d(temp_polyfit_params)
        print(f"Thermodynamic reduced intercept temperature = {temp_polyfit_params[1]} +- {temp_polyfit_errors[1]}")

        """...then plot the data and fitting function"""
        main_legend_1 = axes[1].plot(inverse_squared_log_system_sizes, reduced_intersect_temperatures, marker=".",
                                     markersize=10, color="black", linestyle="None", label="intercept temperatures")
        main_legend_2 = axes[1].plot(continuous_inverse_squared_log_system_sizes,
                                     temp_polynomial_fit(continuous_inverse_squared_log_system_sizes), linestyle="--",
                                     color="black",
                                     label=r"$A + B / (\ln N)^2$ fit $\Rightarrow $ " +
                                           fr"$A =$ {temp_polyfit_params[1]:.3f} $\pm$ {temp_polyfit_errors[1]:.3f}")

        """fit intercept values to power law..."""
        system_sizes = [length ** 2 for index, length in enumerate(linear_system_sizes) if index > 0]
        continuous_system_sizes = np.linspace(0.5 * system_sizes[0], 2.0 * system_sizes[-1], 100)
        value_polyfit_params, value_polyfit_cov = np.polyfit(np.log(system_sizes), np.log(intersect_values), 1,
                                                             cov=True)
        value_polyfit_errors = np.sqrt(np.diag(value_polyfit_cov))
        print(f"Intercept-value power law = C . N^alpha, where "
              f"alpha = {value_polyfit_params[0]} +- {value_polyfit_errors[0]} and "
              f"C = {np.exp(value_polyfit_params[1])} +- {np.exp(value_polyfit_params[1]) * value_polyfit_errors[1]}")

        """...then plot the data and fitting function"""
        inset_legend_1 = inset_axis_2.loglog(system_sizes, intersect_values, marker="*", markersize=10, color="red",
                                             linestyle="None", label="intercept values")
        inset_legend_2 = inset_axis_2.loglog(
            continuous_system_sizes,
            np.exp(value_polyfit_params[1]) * continuous_system_sizes ** value_polyfit_params[0], linestyle="--",
            color="red", label=fr"$C . N^\alpha$ fit $\Rightarrow \alpha =$ {value_polyfit_params[0]:.3f} $\pm$ "
                               fr"{value_polyfit_errors[0]:.3f}")
        supplementary_axis.loglog(system_sizes, intersect_values, marker=".", markersize=10, color="black",
                                  linestyle="None", label="intercept values")
        supplementary_axis.loglog(
            continuous_system_sizes,
            np.exp(value_polyfit_params[1]) * continuous_system_sizes ** value_polyfit_params[0], linestyle="--",
            color="black", label=fr"$C N^\alpha$ fit gives" + "\n" +
                                 fr"$\alpha =$ {value_polyfit_params[0]:.3f} $\pm$ {value_polyfit_errors[0]:.3f} " +
                                 "and" + "\n" + fr"$C =$ {np.exp(value_polyfit_params[1]):.3f} $\pm$ "
                                                fr"{np.exp(value_polyfit_params[1]) * value_polyfit_errors[1]:.3f}")

        """print reduced intercept temperatures and intercept values to file..."""
        intercept_temps_file = open(f"{output_directory}/cramer_von_mises_intercept_temps_and_values_xy_gaussian_noise_"
                                    f"metropolis_w_twists.tsv", "w")
        intercept_temps_file.write("# system size, N".ljust(30) + "reduced intercept temp.".ljust(40) +
                                   "intercept value" + "\n")
        linear_system_sizes.reverse(), linear_system_sizes.pop(), linear_system_sizes.reverse()
        reduced_intersect_temperatures.reverse()
        for system_size_index, length in enumerate(linear_system_sizes):
            intercept_temps_file.write(f"{length}x{length}".ljust(30) +
                                       f"{reduced_intersect_temperatures[system_size_index]:.14e}".ljust(40) +
                                       f"{intersect_values[system_size_index]:.14e}" + "\n")

        """...then print fitting parameters to file"""
        intercept_temps_file.write("\n\n" + f"# Thermodynamic reduced intercept temperature = {temp_polyfit_params[1]} "
                                            f"+- {temp_polyfit_errors[1]}")
        intercept_temps_file.write("\n" +
                                   "# This is a result of fitting A + B / (ln N)^2 to the reduced intercept temps, "
                                   "which gives " + f"A = {temp_polyfit_params[1]} +- {temp_polyfit_errors[1]} and "
                                   f"B = {temp_polyfit_params[0]} +- {temp_polyfit_errors[0]}")
        intercept_temps_file.write(
            "\n\n" + f"# Intercept-value power law = C . N^alpha, where "
                     f"alpha = {value_polyfit_params[0]} +- {value_polyfit_errors[0]} and C = "
                     f"{np.exp(value_polyfit_params[1])} +- {np.exp(value_polyfit_params[1]) * value_polyfit_errors[1]}"
            + "\n" + "# NOTE THAT this power-law fit is OK but not amazing - maybe due to the predicted critical "
                     "slowing down in the phase of the U(1) order parameter")
        intercept_temps_file.close()

        combined_legends_data = main_legend_1 + inset_legend_1 + main_legend_2 + inset_legend_2
        combined_legends_labels = [legend.get_label() for legend in combined_legends_data]
        legends = [axes[0].legend(loc="lower right", ncol=2, fontsize=6.3),
                   axes[1].legend(combined_legends_data, combined_legends_labels, loc="lower right", fontsize=8.75)]
    else:
        legends = [main_axis.legend(loc="lower right", ncol=2, fontsize=11.0)]
        main_axis.text(0.04, 0.37, "dots:  Metrop." + "\n" + "stars: ECMC", fontsize=14.5,
                       transform=main_axis.transAxes,
                       bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(2.4) for legend in legends]
    if symmetry_breaking_paper:
        fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc.pdf",
                    bbox_inches="tight")
    else:
        fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc_"
                    f"sans_intercept_temps.pdf", bbox_inches="tight")

    supplementary_legend = supplementary_axis.legend(loc="lower right", fontsize=12)
    supplementary_legend.get_frame().set_edgecolor("k"), supplementary_legend.get_frame().set_lw(3.0)
    supplementary_fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_"
                              f"intercept_values.pdf", bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


def compute_event_rates(algorithm_name, external_global_moves_string, output_directory, sample_directory, temperatures,
                        length, no_of_samples, no_of_runs, pool):
    try:
        with open(f"{output_directory}/event_rates_{algorithm_name.replace('-', '_')}_{external_global_moves_string}_"
                  f"{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_"
                  f"{no_of_samples}_obs_{no_of_runs}_runs.tsv", "r") as _:
            pass
    except IOError:
        event_rates_file = open(
            f"{output_directory}/event_rates_{algorithm_name.replace('-', '_')}_{external_global_moves_string}_"
            f"{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_{no_of_samples}_"
            f"obs_{no_of_runs}_runs.tsv", "w")
        event_rates_file.write("# temperature".ljust(30) + "event rate".ljust(30) + "event-rate error" + "\n")
        for temperature_index, temperature in enumerate(temperatures):
            event_rate_vs_run = np.array(pool.starmap(get_event_rate,
                                                      [(f"{sample_directory}/run_{run_index}", temperature_index)
                                                       for run_index in range(no_of_runs)]))
            event_rates_file.write(f"{temperature:.14e}".ljust(30) + f"{np.mean(event_rate_vs_run):.14e}".ljust(30) +
                                   f"{np.std(event_rate_vs_run) / len(event_rate_vs_run) ** 0.5:.14e}" + "\n")
        event_rates_file.close()


if __name__ == "__main__":
    if len(sys.argv) > 3:
        raise Exception("InterfaceError: At most two positional arguments are permitted.  None are required but you "
                        "may provide symmetry_breaking_paper (see docstring of main() - default value is True) and "
                        "no_of_system_sizes (which must be an integer greater than 0 and less than 7 - default value "
                        "is 6).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  These must be symmetry_breaking_paper (see docstring of main() - "
              "default value is True) and no_of_system_sizes (which must be an integer greater than 0 and less than 7 "
              "- default value is 6).")
        chosen_no_of_system_sizes = int(sys.argv[2])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 6:
            raise Exception("InterfaceError: chosen_no_of_system_sizes must be an integer greater than 0 and less than "
                            "7 (default value is 6).")
        if sys.argv[1] == "True":
            main(True, chosen_no_of_system_sizes)
        elif sys.argv[1] == "False":
            main(False, chosen_no_of_system_sizes)
        else:
            Exception("InterfaceError: The provided value of symmetry_breaking_paper is neither True nor False (see "
                      "docstring of main()).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be symmetry_breaking_paper (see docstring of main() - "
              "default value is True).  In addition, you may provide number_of_system_sizes - which must be an integer "
              "greater than 0 and less than 7 (default value is 6).")
        if sys.argv[1] == "True":
            main(True)
        elif sys.argv[1] == "False":
            main(False)
        else:
            Exception("InterfaceError: The provided value of symmetry_breaking_paper is neither True nor False (see "
                      "docstring of main()).")
    else:
        print("No positional arguments provided.  None are required but you may provide symmetry_breaking_paper (see "
              "docstring of main() - default value is True) and chosen_no_of_system_sizes (which must be an integer "
              "greater than 0 and less than 7 - default value is 6).")
        main()
