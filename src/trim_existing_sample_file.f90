subroutine trim_existing_sample_file(sample_filename, file_index, temperature_index)
use variables
implicit none
character(100) :: sample_filename, temporary_filename, store_observation
integer :: file_index, temperature_index, observation_index

! if a checkpoint from a previous failed simulation exists, trim the existing sample file to restart from checkpoint

open(file_index, file=sample_filename)
write(temporary_filename, '(A, "/temp_", I2.2, "/temporary.csv")') trim(output_directory), temperature_index
open(file_index + 1, file=temporary_filename)
do observation_index = 0, initial_observation_index
    read(file_index, '(A)') store_observation
    write(file_index + 1, '(A)') store_observation
end do
close(file_index)
close(file_index + 1)
call system('mv ' // trim(temporary_filename) // ' ' // trim(sample_filename))
open(file_index, file=sample_filename, status="old", position="append", action="write")

return
end subroutine trim_existing_sample_file
