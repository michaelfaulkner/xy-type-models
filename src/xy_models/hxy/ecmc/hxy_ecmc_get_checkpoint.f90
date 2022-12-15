subroutine get_checkpoint
use variables
implicit none
integer :: site_index

open(unit=90, file=checkpoint_filename)
read(90, 100) initial_temperature_index
read(90, 100) initial_observation_index
do site_index = 1, no_of_sites
    read(90, 200) spin_field(site_index)
end do
read(90, 100) external_global_moves(1)
read(90, 100) external_global_moves(2)
read(90, 100) no_of_events
read(90, 100) no_of_accepted_external_global_moves
close(90)

call calculate_emergent_field

100 format(I12)
200 format(ES25.14)

return
end subroutine get_checkpoint
