import importlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
polyspectra = importlib.import_module("polyspectra")


def main(config_file, observable_string, no_of_trispectrum_auxiliary_frequency_octaves=4,
         trispectrum_base_period_shift=1):
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
     no_of_jobs, temperature, magnitude_of_temperature_increments, pool) = setup_scripts.set_up_polyspectra_script(
        config_file, observable_string)

    start_time = time.time()
    for _ in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"

        power_trispectrum = polyspectra.get_normalised_power_trispectrum(
            algorithm_name, observable_string, output_directory, temperature_directory, beta, no_of_sites,
            no_of_equilibration_sweeps, no_of_jobs, pool, no_of_trispectrum_auxiliary_frequency_octaves,
            trispectrum_base_period_shift)
        power_trispectrum_zero_mode = polyspectra.get_normalised_power_trispectrum_zero_mode(
            algorithm_name, observable_string, output_directory, temperature_directory, beta, no_of_sites,
            no_of_equilibration_sweeps, no_of_jobs, pool, no_of_trispectrum_auxiliary_frequency_octaves,
            trispectrum_base_period_shift)

        figure, axis = plt.subplots(1, 2, figsize=(10, 5))
        [axis[index].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10) for index in range(2)]
        axis[0].set_ylabel(fr"$|S_X^3 \left( f, f' \right)|$ / $|S_X^3 \left( f_0, f' \right)|$", fontsize=10,
                           labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

        colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_trispectrum_auxiliary_frequency_octaves + 2)))
        for index in range(no_of_trispectrum_auxiliary_frequency_octaves + 2):
            current_color = next(colors)
            if index == 0:
                axis[0].loglog(power_trispectrum_zero_mode[0], power_trispectrum_zero_mode[1], color=current_color,
                               label=r"f' = 0")
            else:
                axis[0].loglog(power_trispectrum[1], power_trispectrum[2][2 ** index - 1], color=current_color,
                               label=fr"f' = {2 ** index}" " x " fr"{power_trispectrum[0][0]:.2e}")

        colors = iter(reversed(plt.cm.rainbow(np.linspace(0, 1, no_of_trispectrum_auxiliary_frequency_octaves + 2))))
        for index in range(no_of_trispectrum_auxiliary_frequency_octaves + 1, -1, -1):
            current_color = next(colors)
            if index == 0:
                axis[1].loglog(power_trispectrum_zero_mode[0], power_trispectrum_zero_mode[1], color=current_color,
                               label=r"f' = 0")
            else:
                axis[1].loglog(power_trispectrum[1], power_trispectrum[2][2 ** index - 1], color=current_color,
                               label=fr"f' = {2 ** index}" " x " fr"{power_trispectrum[0][0]:.2e}")

        figure.tight_layout()
        trispectrum_legend = (axis[0].legend(loc="lower left", fontsize=7.5),
                              axis[1].legend(loc="lower left", fontsize=7.5))
        for legend in trispectrum_legend:
            legend.get_frame().set_edgecolor("k")
            legend.get_frame().set_lw(1.5)
        figure.savefig(f"{output_directory}/{observable_string}_compare_power_trispectrum_auxiliary_frequencies_"
                       f"{no_of_trispectrum_auxiliary_frequency_octaves}_octaves_temp_eq_{temperature:.2f}_"
                       f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                       bbox_inches="tight")
        figure.clf()
        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    if no_of_jobs > 1:
        pool.close()


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate.  In addition, you "
                        "may provide no_of_trispectrum_auxiliary_frequency_octaves (default value is 4) and "
                        "trispectrum_base_period_shift (default value is 1) in the third and fourth positions "
                        "(respectively).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  In addition, you may provide "
              "no_of_trispectrum_auxiliary_frequency_octaves (default value is 4) and trispectrum_base_period_shift "
              "(default value is 1) in the third and fourth positions (respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The third must be no_of_trispectrum_auxiliary_frequency_octaves.  "
              "In addition, you may provide trispectrum_base_period_shift (default value is 1) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The third / fourth must be "
              "no_of_trispectrum_auxiliary_frequency_octaves / trispectrum_base_period_shift.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
