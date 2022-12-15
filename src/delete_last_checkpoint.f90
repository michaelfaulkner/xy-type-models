subroutine delete_last_checkpoint
use variables
implicit none

call system('rm ' // trim(checkpoint_filename))

end subroutine delete_last_checkpoint
