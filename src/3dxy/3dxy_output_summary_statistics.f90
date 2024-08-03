subroutine output_summary_statistics(temperature_index)
use variables
implicit none
character(100) :: temperature_directory, filename
integer :: temperature_index
double precision :: get_monte_carlo_error, magnetic_norm_mean, magnetic_norm_squared_mean, magnetic_norm_quartic_mean
double precision :: magnetic_norm_error, magnetic_susc_mean, magnetic_susc_error
double precision :: potential_mean, potential_squared_mean, potential_error

! make sure the temperature directory (in which to save the summary-stats files) is open
write(temperature_directory, '(A, "/temp_", I2.2)') trim(output_directory), temperature_index
call system('mkdir -p ' // temperature_directory)

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

if (measure_potential) then
    potential_mean = potential_sum / dfloat(no_of_samples)
    potential_squared_mean = potential_squared_sum / dfloat(no_of_samples)
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
