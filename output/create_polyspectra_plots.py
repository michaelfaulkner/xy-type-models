import importlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import time

# import additional modules
setup_scripts = importlib.import_module("setup_scripts")
polyspectra = importlib.import_module("polyspectra")


def main(config_file, observable_string, no_of_trispectrum_auxiliary_frequency_octaves=3,
         trispectrum_base_period_shift=10):
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_temperature_increments,
     no_of_jobs, temperature, magnitude_of_temperature_increments, pool) = setup_scripts.set_up_polyspectra_script(
        config_file, observable_string)

    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))
    figure, axis = plt.subplots(2, figsize=(10, 10))
    axis[1].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10)
    axis[0].set_ylabel(r"$S_X \left( f \right)$ / $S_X \left( f_0 \right)$", fontsize=10, labelpad=10)
    axis[1].set_ylabel(r"$\left| S_X^3 \left( f, f' \right) \right|$ / $\left| S_X^3 \left( f_0, f' \right) \right|$",
                       fontsize=10, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

    start_time = time.time()
    for temperature_index in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")
        current_color = next(colors)
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"

        power_spectrum = polyspectra.get_power_spectrum(algorithm_name, observable_string, output_directory,
                                                        temperature_directory, beta, no_of_sites,
                                                        no_of_equilibration_sweeps, no_of_jobs, pool)
        power_trispectrum = polyspectra.get_power_trispectrum(
            algorithm_name, observable_string, output_directory, temperature_directory, beta, no_of_sites,
            no_of_equilibration_sweeps, no_of_jobs, pool, no_of_trispectrum_auxiliary_frequency_octaves,
            trispectrum_base_period_shift)

        # normalise polyspectra with respect to their low-frequency values
        power_spectrum[1] /= power_spectrum[1, 0]
        power_trispectrum[2] = [spectrum / spectrum[0] for spectrum in power_trispectrum[2]]

        # note frequency indices in order to discard spectrum (trispectrum) for all frequencies >= 0.1 (5.0e-3)  since
        # i) interesting physics at long timescales, and ii) emergent Langevin diffusion breaks down at short timescales
        max_power_spectrum_index = np.argmax(power_spectrum[0] > 1.0e-1) - 1
        max_power_trispectrum_index = np.argmax(power_trispectrum[1] > 5.0e-3) - 1

        axis[0].loglog(power_spectrum[0, :max_power_spectrum_index], power_spectrum[1, :max_power_spectrum_index],
                       color=current_color)
        axis[1].loglog(power_trispectrum[1][:max_power_trispectrum_index],
                       10.0 ** (no_of_temperature_increments - temperature_index) *
                       power_trispectrum[2][len(power_trispectrum[2]) - 1][:max_power_trispectrum_index],
                       color=current_color, label=fr"temperature = {temperature:.2f}; "
                                                  fr"f' = {power_trispectrum[0][len(power_trispectrum[0]) - 1]:.2e}")

        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    figure.tight_layout()
    trispectrum_legend = axis[1].legend(loc="lower left", fontsize=10)
    trispectrum_legend.get_frame().set_edgecolor("k")
    trispectrum_legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/{observable_string}_polyspectra_max_trispectrum_shift_eq_"
                   f"{2 ** no_of_trispectrum_auxiliary_frequency_octaves}_x_{trispectrum_base_period_shift}_delta_t_"
                   f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                   bbox_inches="tight")
    if no_of_jobs > 1:
        pool.close()


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate (in the first and "
                        "second positions, respectively).  In addition, you may provide "
                        "no_of_trispectrum_auxiliary_frequency_octaves (default value is 3) and "
                        "trispectrum_base_period_shift (default value is 10) in the third and fourth positions, "
                        "respectively.")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file / "
              "the string of the observable whose power trispectrum you wish to estimate.  In addition, you may provide"
              " no_of_trispectrum_auxiliary_frequency_octaves (default value is 3) and trispectrum_base_period_shift "
              "(default value is 10) in the third and fourth positions (respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file / the string of the observable whose power trispectrum you wish to estimate / "
              "no_of_trispectrum_auxiliary_frequency_octaves.  In addition, you may provide "
              "trispectrum_base_period_shift (default value is 10) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The first / second / third / fourth must be the location of the "
              "configuration file / the string of the observable whose power trispectrum you wish to estimate / "
              "no_of_trispectrum_auxiliary_frequency_octaves / trispectrum_base_period_shift.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
