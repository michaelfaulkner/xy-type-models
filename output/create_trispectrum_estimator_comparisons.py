import importlib
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path

this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
polyspectra = importlib.import_module("polyspectra")


def main(config_file, observable_string, no_of_trispectrum_octaves=3, trispectrum_base_period_shift=1):
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

        power_trispectrum = polyspectra.try_to_load_normalised_power_trispectrum(beta, no_of_equilibration_sweeps, no_of_jobs,
                                                                     no_of_sites, no_of_trispectrum_octaves,
                                                                     observable_string, output_directory, pool,
                                                                     temperature, temperature_directory,
                                                                     trispectrum_base_period_shift)
        power_trispectrum_direct = polyspectra.try_to_load_normalised_power_trispectrum_direct(beta, no_of_equilibration_sweeps,
                                                                                   no_of_jobs, no_of_sites,
                                                                                   no_of_trispectrum_octaves,
                                                                                   observable_string, output_directory,
                                                                                   pool, temperature,
                                                                                   temperature_directory,
                                                                                   trispectrum_base_period_shift)

        figure, axis = plt.subplots(3, 2, figsize=(10, 10))
        [axis[2, index].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10) for index in range(2)]
        plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

        for index in range(3):
            if index == 0:
                axis[index, 0].loglog(power_trispectrum_direct[1], power_trispectrum_direct[2][index], color='blue',
                                      label="direct estimator")
                axis[index, 0].loglog(power_trispectrum[1], power_trispectrum[2][index], color='red',
                                      label="shortcut estimator")
            elif index == 1:
                axis[index, 0].loglog(power_trispectrum_direct[1], power_trispectrum_direct[2][index], color='blue',
                                      label="direct estimator")
                axis[index, 0].loglog(power_trispectrum[1], power_trispectrum[2][index], color='red',
                                      label="shortcut estimator")
            elif index == 2:
                axis[index, 0].loglog(power_trispectrum_direct[1],
                                      power_trispectrum_direct[2][len(power_trispectrum_direct[2]) - 1],
                                      color='blue', label="direct estimator")
                axis[index, 0].loglog(power_trispectrum[1],
                                      power_trispectrum[2][len(power_trispectrum[2]) - 1], color='red',
                                      label="shortcut estimator")

            if index == 0:
                axis[index, 1].loglog(power_trispectrum[1], power_trispectrum[2][index], color='red',
                                      label="shortcut estimator")
                axis[index, 1].loglog(power_trispectrum_direct[1], power_trispectrum_direct[2][index], color='blue',
                                      label="direct estimator")
            elif index == 1:
                axis[index, 1].loglog(power_trispectrum[1], power_trispectrum[2][index], color='red',
                                      label="shortcut estimator")
                axis[index, 1].loglog(power_trispectrum_direct[1], power_trispectrum_direct[2][index], color='blue',
                                      label="direct estimator")
            elif index == 2:
                axis[index, 1].loglog(power_trispectrum[1],
                                      power_trispectrum[2][len(power_trispectrum[2]) - 1], color='red',
                                      label="shortcut estimator")
                axis[index, 1].loglog(power_trispectrum_direct[1],
                                      power_trispectrum_direct[2][len(power_trispectrum_direct[2]) - 1],
                                      color='blue', label="direct estimator")

        for index in range(3):
            if index == 0:
                axis[index, 0].set_ylabel(fr"$S_X^3 \left( f, f' = 0 \right)$ / $S_X^3 \left( f_0, f'=0 \right)$",
                                          fontsize=10, labelpad=10)
            elif index == 1:
                axis[index, 0].set_ylabel(fr"$|S_X^3 \left( f, f_0' \right)|$ / $|S_X^3 \left( f_0, f_0' \right)|$, "
                                          fr"$f_0' = {power_trispectrum[0][0]:.2e}$", fontsize=10, labelpad=10)
            else:
                axis[index, 0].set_ylabel(fr"$|S_X^3 \left( f, {2 ** no_of_trispectrum_octaves} f_0' \right)|$ / "
                                          fr"$|S_X^3 \left(f_0, {2 ** no_of_trispectrum_octaves} f_0' \right)|$, "
                                          fr"$f_0' = {power_trispectrum[0][0]:.2e}$", fontsize=10, labelpad=10)

        figure.tight_layout()
        trispectrum_legend = (axis[0, 0].legend(loc="lower left", fontsize=10),
                              axis[0, 1].legend(loc="lower left", fontsize=10))
        for legend in trispectrum_legend:
            legend.get_frame().set_edgecolor("k")
            legend.get_frame().set_lw(1.5)
        figure.savefig(f"{output_directory}/{observable_string}_compare_power_trispectrum_functions_"
                       f"{no_of_trispectrum_octaves}_octaves_temp_eq_{temperature:.2f}.pdf", bbox_inches="tight")
        figure.clf()
        temperature -= magnitude_of_temperature_increments

    if no_of_jobs > 1:
        pool.close()


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate.  In addition, you "
                        "may provide no_of_trispectrum_octaves (default value is 3) and trispectrum_base_period_shift "
                        "(default value is 1).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  In addition, you may provide no_of_trispectrum_octaves (default "
              "value is 3) and trispectrum_base_period_shift (default value is 1).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The third must be no_of_trispectrum_octaves.  In addition, you may"
              " provide trispectrum_base_period_shift (default value is 1).")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The third / fourth must be no_of_trispectrum_octaves / "
              "trispectrum_base_period_shift.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
