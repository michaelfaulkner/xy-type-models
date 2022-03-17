subroutine create_sample_file(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index, file_index

write(filename, '(A, "/temp_", I2.2, "_sample_potential.csv")') trim(output_directory), temperature_index
open(unit=20, file=filename)
if (use_external_global_moves) then
    write(filename, '(A, "/temp_", I2.2, "_sample_external_global_move.csv")') trim(output_directory), temperature_index
    open(unit=30, file=filename)
end if
write(filename, '(A, "/temp_", I2.2, "_sample_non_normalised_magnetisation.csv")') trim(output_directory), temperature_index
open(unit=40, file=filename)
write(filename, '(A, "/temp_", I2.2, "_sample_potential_1st_deriv_of_potential.csv")') trim(output_directory), temperature_index
open(unit=50, file=filename)
write(filename, '(A, "/temp_", I2.2, "_sample_potential_2nd_deriv_of_potential.csv")') trim(output_directory), temperature_index
open(unit=60, file=filename)
if (calculate_potential_minimising_twists) then
    write(filename, '(A, "/temp_", I2.2, "_sample_potential_minimising_twists.csv")') trim(output_directory), temperature_index
    open(unit=70, file=filename)
end if

do file_index = 2, 7
    write(10 * file_index, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
                                    trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature
end do

return
end subroutine create_sample_file
