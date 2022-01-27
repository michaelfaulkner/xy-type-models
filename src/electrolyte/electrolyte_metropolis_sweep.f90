subroutine metropolis_sweep
use variables
implicit none
integer :: i, lattice_site
double precision :: attempt_charge_hop_or_field_rotation

call randomise_array_of_sites
do i = 1, no_of_sites
    lattice_site = array_of_sites(i)
    attempt_charge_hop_or_field_rotation = rand()
    if (attempt_charge_hop_or_field_rotation < charge_hop_proportion_over_two) then
        call attempt_charge_hop(lattice_site, 1)
    else if (attempt_charge_hop_or_field_rotation < charge_hop_proportion) then
        call attempt_charge_hop(lattice_site, 2)
    else
        call attempt_field_rotation(lattice_site)
    end if
end do

return
end subroutine metropolis_sweep


subroutine attempt_field_rotation(lattice_site)
use variables
implicit none
integer :: lattice_site
double precision :: field_rotation_value, potential_difference, candidate_electric_field_component(4)

field_rotation_value = width_of_proposal_interval * (rand() - 0.5d0)
candidate_electric_field_component(1) = electric_field(lattice_site, 1) + field_rotation_value
candidate_electric_field_component(2) = electric_field(lattice_site, 2) - field_rotation_value
candidate_electric_field_component(3) = electric_field(get_north_neighbour(lattice_site), 1) - field_rotation_value
candidate_electric_field_component(4) = electric_field(get_east_neighbour(lattice_site), 2) + field_rotation_value

potential_difference = 0.5d0 &
                        * (candidate_electric_field_component(1) ** 2 + candidate_electric_field_component(2) ** 2 &
                         + candidate_electric_field_component(3) ** 2 + candidate_electric_field_component(4) ** 2 &
                         - electric_field(lattice_site, 1) ** 2 - electric_field(lattice_site, 2) ** 2 &
                         - electric_field(get_north_neighbour(lattice_site), 1) ** 2 &
                         - electric_field(get_east_neighbour(lattice_site), 2) ** 2)

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    electric_field(lattice_site, 1) = candidate_electric_field_component(1)
    electric_field(lattice_site, 2) = candidate_electric_field_component(2)
    electric_field(get_north_neighbour(lattice_site), 1) = candidate_electric_field_component(3)
    electric_field(get_east_neighbour(lattice_site), 2) = candidate_electric_field_component(4)
    no_of_accepted_field_rotations = no_of_accepted_field_rotations + 1
end if
  
return
end subroutine attempt_field_rotation
