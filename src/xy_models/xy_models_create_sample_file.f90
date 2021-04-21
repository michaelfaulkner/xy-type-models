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
open(unit=10, file=filename)

write(10, '(A1, A23, "; ", I0.4, " x ", I0.4, " lattice sites; temperature = ", ES8.2)') "#", algorithm_name, side, &
                                                                                                side, temperature
write(10, '(A1, A29, 8A30)') "#", "potential", "magnetisation_x", "magnetisation_y", &
                                "1st_deriv_of_potential_x", "1st_deriv_of_potential_y", &
                                "2nd_deriv_of_potential_x", "2nd_deriv_of_potential_y", &
                                "minimising_twist_field_x", "minimising_twist_field_y"

return
end subroutine create_sample_file