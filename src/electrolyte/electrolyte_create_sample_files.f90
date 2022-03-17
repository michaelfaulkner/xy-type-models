subroutine create_sample_files(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index, file_index

write(filename, '(A, "/temp_", I2.2, "_potential.csv")') trim(output_directory), temperature_index
open(unit=20, file=filename)
write(filename, '(A, "/temp_", I2.2, "_electric_field_sum.csv")') trim(output_directory), temperature_index
open(unit=30, file=filename)
write(filename, '(A, "/temp_", I2.2, "_external_global_moves.csv")') trim(output_directory), temperature_index
open(unit=40, file=filename)

do file_index = 2, 4
    write(10 * file_index, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
                                    trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature
end do

return
end subroutine create_sample_files
