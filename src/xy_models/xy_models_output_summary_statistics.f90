subroutine output_summary_statistics(temperature_index)
use variables
implicit none
character(100) ::  filename
integer :: temperature_index
double precision :: get_monte_carlo_error, magnetic_norm_mean, magnetic_norm_squared_mean, magnetic_norm_quartic_mean
double precision :: magnetic_norm_error, magnetic_susc_mean, magnetic_susc_error, inverse_vacuum_perm_mean
double precision :: inverse_vacuum_perm_squared_mean, inverse_vacuum_perm_error, macro_josephson_current_mean(2)
double precision :: macro_josephson_current_squared_mean, macro_josephson_current_quartic_mean
double precision :: macro_josephson_current_susc_mean, macro_josephson_current_susc_error, helicity_mean, helicity_error
double precision :: potential_mean, potential_squared_mean, potential_quartic_mean, potential_error
double precision :: specific_heat_per_particle_mean, specific_heat_per_particle_error, hot_twist_relaxations_mean(2)
double precision :: hot_twist_relaxations_squared_mean, hot_twist_relaxations_quartic_mean
double precision :: potential_minimising_twist_susc_mean, potential_minimising_twist_susc_error
double precision :: twist_relaxations_mean(2), twist_relaxations_squared_mean, twist_relaxations_quartic_mean
double precision :: twist_relaxation_susc_mean, twist_relaxation_susc_error, zero_mode_mean(2), zero_mode_squared_mean
double precision :: zero_mode_quartic_mean, zero_mode_susc_mean, zero_mode_susc_error, topological_sector_mean(2)
double precision :: topological_sector_squared_mean, topological_sector_quartic_mean, topological_susc_mean
double precision :: topological_susc_error

if (measure_magnetisation) then
    ! magnetic_norm_squared = raw_magnetic_norm_squared / no_of_sites ** 2
    magnetic_norm_mean = raw_magnetic_norm_sum / dfloat(no_of_sites * no_of_samples)
    magnetic_norm_squared_mean = raw_magnetic_norm_squared_sum / dfloat(no_of_sites ** 2 * no_of_samples)
    magnetic_norm_quartic_mean = raw_magnetic_norm_quartic_sum / dfloat(no_of_sites ** 4 * no_of_samples)
    magnetic_norm_error = get_monte_carlo_error(magnetic_norm_mean, magnetic_norm_squared_mean)
    magnetic_susc_mean = dfloat(no_of_sites) * beta * (magnetic_norm_squared_mean - magnetic_norm_mean ** 2)
    magnetic_susc_error = dfloat(no_of_sites) * beta * &
                            get_monte_carlo_error(magnetic_norm_squared_mean, magnetic_norm_quartic_mean)

    write(filename, '(A, "/temp_", I2.2, "/magnetic_summary_stats.csv")') trim(output_directory), temperature_index
    open(unit=11, file=filename)
    write(11, 200) magnetic_norm_mean
    write(11, 200) magnetic_norm_error
    write(11, 200) magnetic_susc_mean
    write(11, 200) magnetic_susc_error
    close(11)
end if

if (measure_helicity) then
    ! inverse_vacuum_perm = raw_inverse_vacuum_perm / 2 / no_of_sites
    inverse_vacuum_perm_mean = 0.5d0 * raw_inverse_vacuum_perm_sum / dfloat(no_of_sites * no_of_samples)
    inverse_vacuum_perm_squared_mean = 0.25d0 * raw_inverse_vacuum_perm_squared_sum / &
                                                dfloat(no_of_sites ** 2 * no_of_samples)
    inverse_vacuum_perm_error = get_monte_carlo_error(inverse_vacuum_perm_mean, inverse_vacuum_perm_squared_mean)

    ! macro_josephson_current = raw_macro_josephson_current / no_of_sites
    macro_josephson_current_mean(1) = raw_macro_josephson_current_sum(1) / dfloat(no_of_sites * no_of_samples)
    macro_josephson_current_mean(2) = raw_macro_josephson_current_sum(2) / dfloat(no_of_sites * no_of_samples)
    macro_josephson_current_squared_mean = raw_macro_josephson_current_squared_sum / dfloat(no_of_sites ** 2 * no_of_samples)
    macro_josephson_current_quartic_mean = raw_macro_josephson_current_quartic_sum / dfloat(no_of_sites ** 4 * no_of_samples)
    macro_josephson_current_susc_mean = no_of_sites * beta * 0.5d0 * (macro_josephson_current_squared_mean - &
                                            macro_josephson_current_mean(1) ** 2 - macro_josephson_current_mean(2) ** 2)
    macro_josephson_current_susc_error = no_of_sites * beta * 0.5d0 * &
                    get_monte_carlo_error(macro_josephson_current_squared_mean, macro_josephson_current_quartic_mean)

    helicity_mean = inverse_vacuum_perm_mean - macro_josephson_current_susc_mean
    helicity_error = (inverse_vacuum_perm_error ** 2 + macro_josephson_current_susc_error ** 2) ** 0.5

    write(filename, '(A, "/temp_", I2.2, "/helicity_summary_stats.csv")') trim(output_directory), temperature_index
    open(unit=11, file=filename)
    write(11, 200) helicity_mean
    write(11, 200) helicity_error
    close(11)
end if

if (measure_potential) then
    potential_mean = potential_sum / dfloat(no_of_samples)
    potential_squared_mean = potential_squared_sum / dfloat(no_of_samples)
    potential_quartic_mean = potential_quartic_sum / dfloat(no_of_samples)
    potential_error = get_monte_carlo_error(potential_mean, potential_squared_mean)
    specific_heat_per_particle_mean = beta ** 2 * (potential_squared_mean - potential_mean ** 2) / dfloat(no_of_sites)
    specific_heat_per_particle_error = beta ** 2 * &
                            get_monte_carlo_error(potential_squared_mean, potential_quartic_mean) / dfloat(no_of_sites)

    write(filename, '(A, "/temp_", I2.2, "/potential_summary_stats.csv")') trim(output_directory), temperature_index
    open(unit=11, file=filename)
    write(11, 200) potential_mean
    write(11, 200) potential_error
    write(11, 200) specific_heat_per_particle_mean
    write(11, 200) specific_heat_per_particle_error
    close(11)
end if

if (measure_hot_twist_relaxations) then
    hot_twist_relaxations_mean(1) = dfloat(hot_twist_relaxations_sum(1)) / dfloat(no_of_samples)
    hot_twist_relaxations_mean(2) = dfloat(hot_twist_relaxations_sum(2)) / dfloat(no_of_samples)
    hot_twist_relaxations_squared_mean = dfloat(hot_twist_relaxations_squared_sum) / dfloat(no_of_samples)
    hot_twist_relaxations_quartic_mean = dfloat(hot_twist_relaxations_quartic_sum) / dfloat(no_of_samples)

    potential_minimising_twist_susc_mean = 2.0d0 * pi ** 2 * beta * (hot_twist_relaxations_squared_mean - &
                                                hot_twist_relaxations_mean(1) ** 2 - hot_twist_relaxations_mean(2) ** 2)
    potential_minimising_twist_susc_error = 2.0d0 * pi ** 2 * beta * &
                            get_monte_carlo_error(hot_twist_relaxations_squared_mean, hot_twist_relaxations_quartic_mean)

    write(filename, '(A, "/temp_", I2.2, "/hot_twist_relaxations_summary_stats.csv")') trim(output_directory), &
                                                                                            temperature_index
    open(unit=11, file=filename)
    write(11, 200) potential_minimising_twist_susc_mean
    write(11, 200) potential_minimising_twist_susc_error
    close(11)
end if

if (measure_twist_relaxations) then
    twist_relaxations_mean(1) = dfloat(twist_relaxations_sum(1)) / dfloat(no_of_samples)
    twist_relaxations_mean(2) = dfloat(twist_relaxations_sum(2)) / dfloat(no_of_samples)
    twist_relaxations_squared_mean = dfloat(twist_relaxations_squared_sum) / dfloat(no_of_samples)
    twist_relaxations_quartic_mean = dfloat(twist_relaxations_quartic_sum) / dfloat(no_of_samples)

    twist_relaxation_susc_mean = 2.0d0 * pi ** 2 * beta * &
                    (twist_relaxations_squared_mean - twist_relaxations_mean(1) ** 2 - twist_relaxations_mean(2) ** 2)
    twist_relaxation_susc_error = 2.0d0 * pi ** 2 * beta * &
                    get_monte_carlo_error(twist_relaxations_squared_mean, twist_relaxations_quartic_mean)

    write(filename, '(A, "/temp_", I2.2, "/twist_relaxations_summary_stats.csv")') trim(output_directory), &
                                                                                        temperature_index
    open(unit=11, file=filename)
    write(11, 200) twist_relaxation_susc_mean
    write(11, 200) twist_relaxation_susc_error
    close(11)
end if

if (measure_emergent_field) then
    zero_mode_mean(1) = sum_of_emergent_field_sum(1) / dfloat(no_of_sites * no_of_samples)
    zero_mode_mean(2) = sum_of_emergent_field_sum(2) / dfloat(no_of_sites * no_of_samples)
    zero_mode_squared_mean = sum_of_emergent_field_squared_sum / dfloat(no_of_sites ** 2 * no_of_samples)
    zero_mode_quartic_mean = sum_of_emergent_field_quartic_sum / dfloat(no_of_sites ** 4 * no_of_samples)

    zero_mode_susc_mean = 0.5d0 * dfloat(no_of_sites) * beta * &
                                (zero_mode_squared_mean - zero_mode_mean(1) ** 2 - zero_mode_mean(2) ** 2)
    zero_mode_susc_error = 0.5d0 * dfloat(no_of_sites) * beta * &
                                get_monte_carlo_error(zero_mode_squared_mean, zero_mode_quartic_mean)

    topological_sector_mean(1) = dfloat(topological_sector_sum(1)) / dfloat(no_of_samples)
    topological_sector_mean(2) = dfloat(topological_sector_sum(2)) / dfloat(no_of_samples)
    topological_sector_squared_mean = dfloat(topological_sector_squared_sum) / dfloat(no_of_samples)
    topological_sector_quartic_mean = dfloat(topological_sector_quartic_sum) / dfloat(no_of_samples)

    topological_susc_mean = 2.0d0 * pi ** 2 * beta * &
                (topological_sector_squared_mean - topological_sector_mean(1) ** 2 - topological_sector_mean(2) ** 2)
    topological_susc_error = 2.0d0 * pi ** 2 * beta * &
                                get_monte_carlo_error(topological_sector_squared_mean, topological_sector_quartic_mean)

    write(filename, '(A, "/temp_", I2.2, "/emergent_field_summary_stats.csv")') trim(output_directory), temperature_index
    open(unit=11, file=filename)
    write(11, 200) zero_mode_susc_mean
    write(11, 200) zero_mode_susc_error
    write(11, 200) topological_susc_mean
    write(11, 200) topological_susc_error
    close(11)
end if

200 format(ES25.14)

return
end subroutine output_summary_statistics
