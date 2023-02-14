subroutine get_checkpoint
use variables
implicit none
integer :: site_index

open(unit=11, file=checkpoint_filename)
read(11, 100) initial_temperature_index
read(11, 100) initial_observation_index
read(11, 200) width_of_proposal_interval
do site_index = 1, no_of_sites
    read(11, 200) electric_field(site_index, 1)
    read(11, 200) electric_field(site_index, 2)
end do
read(11, 100) net_charge_displacement(1)
read(11, 100) net_charge_displacement(2)
read(11, 100) external_global_moves(1)
read(11, 100) external_global_moves(2)
read(11, 200) no_of_accepted_charge_hops_per_site
read(11, 200) no_of_accepted_field_rotations_per_site
read(11, 200) no_of_accepted_external_global_moves
close(11)

100 format(I12)
200 format(ES25.14)

return
end subroutine get_checkpoint
