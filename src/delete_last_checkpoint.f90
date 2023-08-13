subroutine delete_last_checkpoint
use variables
implicit none
character(100) :: temporary_filename
logical :: checkpoint_exists

! create end_of_simulation.csv to skip this simulation when restarting multiple parallel simulations that had been -
! simulations restart from the checkpoint of some run of larger index that has not finished
write(temporary_filename, '(A, "/end_of_simulation.csv")') trim(output_directory)
open(unit=11, file=temporary_filename)
write(11, '(A)') 'simulation is complete'
close(11)

write(checkpoint_filename, '(A, "/checkpoint.csv")') trim(output_directory)
inquire(file=checkpoint_filename, exist=checkpoint_exists)
if (checkpoint_exists) then
    call system('rm ' // trim(checkpoint_filename))
end if

end subroutine delete_last_checkpoint
