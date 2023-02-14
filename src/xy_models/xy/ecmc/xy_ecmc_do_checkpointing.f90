subroutine do_checkpointing(temperature_index, observation_index)
use variables
implicit none
character(100) :: temporary_filename
integer :: temperature_index, observation_index, site_index
real :: current_time

! a checkpoint is a snapshot of the system, created in case the simulation fails - this subroutine creates a new one if
! the elapsed CPU time since the last checkpoint is greater than time_between_checkpoints

write(temporary_filename, '(A, "/temporary_checkpoint.csv")') trim(output_directory)
call cpu_time(current_time)
if (current_time - previous_checkpointing_time > time_between_checkpoints) then
    open(unit=11, file=temporary_filename)
    write(11, 100) temperature_index
    write(11, 100) observation_index
    do site_index = 1, no_of_sites
        write(11, 200) spin_field(site_index)
    end do
    write(11, 100) external_global_moves(1)
    write(11, 100) external_global_moves(2)
    write(11, 200) no_of_events_per_unit_spin_space_distance
    write(11, 200) no_of_accepted_external_global_moves
    close(11)

    call system('mv ' // trim(temporary_filename) // ' ' // trim(checkpoint_filename))
    previous_checkpointing_time = current_time
end if

100 format(I12)
200 format(ES25.14)

return
end subroutine do_checkpointing
