subroutine create_sample_files
use variables
implicit none
character(100) :: filename
character(100) :: temperature_directory
character(100), parameter :: temperature_string = "/temp_eq_"

! OPENS NEW DIRECTORY IN WHICH TO SAVE THE MARKOV CHAIN FOR THE CURRENT TEMPERATURE
write (temperature_directory, '(A, F4.2)' ) trim(output_directory)//trim(temperature_string), temperature
call system ( 'mkdir -p ' // temperature_directory )

write (filename, '(A, F4.2, "//magnetisation_x_sample.dat")' ) trim(output_directory)//trim(temperature_string), &
        temperature
open(unit=10,file=filename)
write (filename, '(A, F4.2, "//magnetisation_y_sample.dat")' ) trim(output_directory)//trim(temperature_string), &
        temperature
open(unit=11,file=filename)
write (filename, '(A, F4.2, "//mean_1st_derivative_of_potential_x_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), temperature
open(unit=12,file=filename)
write (filename, '(A, F4.2, "//mean_1st_derivative_of_potential_y_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), temperature
open(unit=13,file=filename)
write (filename, '(A, F4.2, "//mean_2nd_derivative_of_potential_x_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), temperature
open(unit=14,file=filename)
write (filename, '(A, F4.2, "//mean_2nd_derivative_of_potential_y_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), temperature
open(unit=15,file=filename)
write (filename, '(A, F4.2, "//potential_sample.dat")' ) trim(output_directory)//trim(temperature_string), temperature
open(unit=16, file = filename)
write (filename, '(A, F4.2, "//external_minimising_twist_field_x_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), temperature
open(unit=17, file = filename)
write (filename, '(A, F4.2, "//external_minimising_twist_field_y_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), temperature
open(unit=18, file = filename)

return
end subroutine create_sample_files
