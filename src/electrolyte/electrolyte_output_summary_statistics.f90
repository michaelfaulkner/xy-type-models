subroutine output_summary_statistics(temperature_index)
use variables
implicit none
character(100) :: temperature_directory, filename
integer :: temperature_index
double precision :: get_monte_carlo_error, potential_mean, potential_squared_mean, potential_error, zero_mode_mean(2)
double precision :: zero_mode_squared_mean, zero_mode_quartic_mean, zero_mode_susc_mean, zero_mode_susc_error
double precision :: topological_sector_mean(2), topological_sector_squared_mean, topological_sector_quartic_mean
double precision :: topological_susc_mean, topological_susc_error

! make sure the temperature directory (in which to save the summary-stats files) is open
write(temperature_directory, '(A, "/temp_", I2.2)') trim(output_directory), temperature_index
call system('mkdir -p ' // temperature_directory)

if (measure_electric_field_sum) then
    ! raw_electric_field_zero_mode = \sum_r E(r) / two_pi = no_of_sites * Ebar / c
    zero_mode_mean(1) =  two_pi * dfloat(raw_electric_field_zero_mode_sum(1)) / dfloat(no_of_sites * no_of_observations)
    zero_mode_mean(2) = two_pi * dfloat(raw_electric_field_zero_mode_sum(2)) / dfloat(no_of_sites * no_of_observations)
    zero_mode_squared_mean = two_pi ** 2 * dfloat(raw_electric_field_zero_mode_squared_sum) / &
                                                                        dfloat(no_of_sites ** 2 * no_of_observations)
    zero_mode_quartic_mean = two_pi ** 4 * dfloat(raw_electric_field_zero_mode_quartic_sum) / &
                                                                        dfloat(no_of_sites ** 4 * no_of_observations)

    zero_mode_susc_mean = 0.5d0 * no_of_sites * beta * &
                                (zero_mode_squared_mean - zero_mode_mean(1) ** 2 - zero_mode_mean(2) ** 2)
    zero_mode_susc_error = 0.5d0 * no_of_sites * beta * &
                                get_monte_carlo_error(zero_mode_squared_mean, zero_mode_quartic_mean)

    topological_sector_mean(1) = dfloat(topological_sector_sum(1)) / dfloat(no_of_observations)
    topological_sector_mean(2) = dfloat(topological_sector_sum(2)) / dfloat(no_of_observations)
    topological_sector_squared_mean = dfloat(topological_sector_squared_sum) / dfloat(no_of_observations)
    topological_sector_quartic_mean = dfloat(topological_sector_quartic_sum) / dfloat(no_of_observations)

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
    potential_mean = potential_sum / dfloat(no_of_observations)
    potential_squared_mean = potential_squared_sum / dfloat(no_of_observations)
    potential_error = get_monte_carlo_error(potential_mean, potential_squared_mean)

    write(filename, '(A, "/temp_", I2.2, "/potential_summary_stats.csv")') trim(output_directory), temperature_index
    open(unit=11, file=filename)
    write(11, 200) potential_mean
    write(11, 200) potential_error
    close(11)
end if

200 format(ES25.14)

return
end subroutine output_summary_statistics
