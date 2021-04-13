subroutine create_sample_files
use variables
implicit none
character(100) :: filename
character(100) :: temperature_directory
character(100), parameter :: temperature_string = "/temp_eq_"

! OPENS NEW DIRECTORY IN WHICH TO SAVE THE MARKOV CHAIN FOR THE CURRENT TEMPERATURE
write (temperature_directory, '(A, F4.2)' ) trim(output_directory)//trim(temperature_string), temperature
call system ( 'mkdir -p ' // temperature_directory )

write(filename, '(A, F4.2, "//field_sum_sample.dat")') trim(output_directory)//trim(temperature_string), temperature
open(unit=10, file=filename)
write(filename, '(A, F4.2, "//potential_sample.dat")') trim(output_directory)//trim(temperature_string), temperature
open(unit=11, file=filename)

return
end subroutine create_sample_files
