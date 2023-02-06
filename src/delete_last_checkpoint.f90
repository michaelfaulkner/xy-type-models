subroutine delete_last_checkpoint
use variables
implicit none
character(100) :: temporary_filename

write(temporary_filename, '(A, "/end_of_simulation.csv")') trim(output_directory)
open(unit=90, file=temporary_filename)
write(90, '(A)') 'simulation is complete'
close(90)
call system('rm ' // trim(checkpoint_filename))

end subroutine delete_last_checkpoint
