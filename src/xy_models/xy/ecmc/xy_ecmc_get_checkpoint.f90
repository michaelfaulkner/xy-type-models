subroutine get_checkpoint
use variables
implicit none
integer :: site_index

open(unit=11, file=checkpoint_filename)
read(11, 100) initial_temperature_index
read(11, 100) initial_sample_index
do site_index = 1, no_of_sites
    read(11, 200) spin_field(site_index)
end do
read(11, 100) external_global_moves(1)
read(11, 100) external_global_moves(2)
read(11, 200) no_of_events_per_unit_spin_space_distance
read(11, 200) no_of_accepted_external_global_moves
close(11)

100 format(I12)
200 format(ES25.14)

return
end subroutine get_checkpoint
