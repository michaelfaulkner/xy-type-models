subroutine create_sample_file
use variables
implicit none
character(100) :: filename
character(100) :: temperature_directory
character(100), parameter :: temperature_string = "/temp_eq_"

! opens new directory in which to save the sample files at the current temperature
write(temperature_directory, '(A, F4.2)') trim(output_directory)//trim(temperature_string), temperature
call system('mkdir -p ' // temperature_directory)

write(filename, '(A, F4.2, "//sample.csv")') trim(output_directory)//trim(temperature_string), temperature
open(unit=20, file=filename)

write(20, '(A1, A23, "; ", I0.4, " x ", I0.4, " lattice sites; temperature = ", ES8.2)') "#", algorithm_name, side, &
                                                                                                side, temperature
write(20, '(A1, A29, 2A30)') "#", "potential", "sum_electric_field_x", "sum_electric_field_y"

return
end subroutine create_sample_file
