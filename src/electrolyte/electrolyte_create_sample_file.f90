subroutine create_sample_file(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index

write(filename, '(A, "/temp_", I2.2, "_sample.csv")') trim(output_directory), temperature_index
open(unit=20, file=filename)

write(20, '(A1, A23, "; ", I0.4, " x ", I0.4, " lattice no_of_sites; temperature = ", ES8.2)') "#", &
                                        algorithm_name, integer_lattice_length, integer_lattice_length, temperature
write(20, '(A1, A29, 4A30)') "#", "potential", "external_global_move_x", "external_global_move_y", &
                                    "sum_of_electric_field_x", "sum_of_electric_field_y"


return
end subroutine create_sample_file
