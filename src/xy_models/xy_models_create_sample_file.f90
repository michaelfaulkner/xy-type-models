subroutine create_sample_file(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index

write(filename, '(A, "/temp_", I2.2, "_sample.csv")') trim(output_directory), temperature_index
open(unit=20, file=filename)

write(20, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') trim(algorithm_name), &
                integer_lattice_length, integer_lattice_length, temperature
write(20, '("#", A29, 10A30)') "potential", "external_global_move_x", "external_global_move_y", &
                                "nonnormalised magnetisation_x", "nonnormalised magnetisation_y", &
                                "1st_deriv_of_potential_x", "1st_deriv_of_potential_y", &
                                "2nd_deriv_of_potential_x", "2nd_deriv_of_potential_y", &
                                "potential_minimising_twists_x", "potential_minimising_twists_y"

return
end subroutine create_sample_file
