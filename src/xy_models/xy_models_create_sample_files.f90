subroutine create_sample_files(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index, file_index

write(filename, '(A, "/temp_", I2.2, "_potential.csv")') trim(output_directory), temperature_index
open(unit=20, file=filename)
write(filename, '(A, "/temp_", I2.2, "_magnetisation.csv")') trim(output_directory), temperature_index
open(unit=30, file=filename)
write(filename, '(A, "/temp_", I2.2, "_1st_deriv_of_potential.csv")') trim(output_directory), temperature_index
open(unit=40, file=filename)
write(filename, '(A, "/temp_", I2.2, "_2nd_deriv_of_potential.csv")') trim(output_directory), temperature_index
open(unit=50, file=filename)

do file_index = 2, 5
    write(10 * file_index, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
                                    trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature
end do

if (use_external_global_moves) then
    write(filename, '(A, "/temp_", I2.2, "_external_global_moves.csv")') trim(output_directory), temperature_index
    open(unit=60, file=filename)
    write(60, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
                                    trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature
end if
if (calculate_potential_minimising_twists) then
    write(filename, '(A, "/temp_", I2.2, "_potential_minimising_twists.csv")') trim(output_directory), temperature_index
    open(unit=70, file=filename)
    write(70, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
                                    trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature
end if

return
end subroutine create_sample_files
