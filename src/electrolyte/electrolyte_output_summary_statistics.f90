subroutine output_summary_statistics(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index
double precision :: get_monte_carlo_error, potential_mean, potential_squared_mean, potential_quartic_mean, potential_error
double precision :: specific_heat_per_particle_mean, specific_heat_per_particle_error, zero_mode_mean(2)
double precision :: zero_mode_squared_mean, zero_mode_quartic_mean, zero_mode_susc_mean, zero_mode_susc_error
double precision :: topological_sector_mean(2), topological_sector_squared_mean, topological_sector_quartic_mean
double precision :: topological_susc_mean, topological_susc_error

if (measure_electric_field_sum) then
    ! raw_electric_field_zero_mode = \sum_r E(r) / no_of_sites / two_pi = Ebar / two_pi
    zero_mode_mean(1) =  two_pi * raw_electric_field_zero_mode_sum(1) / dfloat(no_of_samples)
    zero_mode_mean(2) = two_pi * raw_electric_field_zero_mode_sum(2) / dfloat(no_of_samples)
    zero_mode_squared_mean = two_pi ** 2 * raw_electric_field_zero_mode_squared_sum / dfloat(no_of_samples)
    zero_mode_quartic_mean = two_pi ** 4 * raw_electric_field_zero_mode_quartic_sum / dfloat(no_of_samples)

    zero_mode_susc_mean = 0.5d0 * dfloat(no_of_sites) * beta * &
                                (zero_mode_squared_mean - zero_mode_mean(1) ** 2 - zero_mode_mean(2) ** 2)
    zero_mode_susc_error = 0.5d0 * dfloat(no_of_sites) * beta * &
                                get_monte_carlo_error(zero_mode_squared_mean, zero_mode_quartic_mean)

    topological_sector_mean(1) = dfloat(topological_sector_sum(1)) / dfloat(no_of_samples)
    topological_sector_mean(2) = dfloat(topological_sector_sum(2)) / dfloat(no_of_samples)
    topological_sector_squared_mean = dfloat(topological_sector_squared_sum) / dfloat(no_of_samples)
    topological_sector_quartic_mean = dfloat(topological_sector_quartic_sum) / dfloat(no_of_samples)

    topological_susc_mean = 0.5d0 * two_pi ** 2 * beta * &
                (topological_sector_squared_mean - topological_sector_mean(1) ** 2 - topological_sector_mean(2) ** 2)
    topological_susc_error = 0.5d0 * two_pi ** 2 * beta * &
                                get_monte_carlo_error(topological_sector_squared_mean, topological_sector_quartic_mean)

    write(filename, '(A, "/temp_", I2.2, "/electric_field_summary_stats.csv")') trim(output_directory), temperature_index
    open(unit=11, file=filename)
    write(11, 200) zero_mode_susc_mean
    write(11, 200) zero_mode_susc_error
    write(11, 200) topological_susc_mean
    write(11, 200) topological_susc_error
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

200 format(ES25.14)

return
end subroutine output_summary_statistics
