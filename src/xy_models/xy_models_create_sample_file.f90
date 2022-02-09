subroutine create_sample_file
use variables
implicit none
character(100) :: filename
character(100) :: temperature_directory
character(100), parameter :: temperature_string = "/temp_eq_"

! opens new directory in which to save the sample files at the current temperature
write(temperature_directory, '(A, F4.2)') trim(output_directory)//trim(temperature_string), temperature
call system('mkdir -p ' // temperature_directory)

write(filename, '(A, F4.2, "//sample.csv")') trim(output_directory)//trim(temperature_string), &
        temperature
open(unit=20, file=filename)

write(20, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') trim(algorithm_name), &
                integer_lattice_length, integer_lattice_length, temperature
write(20, '("#", A29, 10A30)') "potential", "external_global_move_x", "external_global_move_y", &
                                "nonnormalised magnetisation_x", "nonnormalised magnetisation_y", &
                                "1st_deriv_of_potential_x", "1st_deriv_of_potential_y", &
                                "2nd_deriv_of_potential_x", "2nd_deriv_of_potential_y", &
                                "minimising_twist_field_x", "minimising_twist_field_y"

return
end subroutine create_sample_file
