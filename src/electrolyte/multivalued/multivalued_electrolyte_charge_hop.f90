subroutine charge_hop(lattice_site, cartesian_component)
use variables
implicit none
integer :: lattice_site, cartesian_component
double precision :: potential_difference, candidate_electric_field_component

candidate_electric_field_component = electric_field(lattice_site, cartesian_component) + elementary_charge &
                                        * (2.0d0 * floor(2.0d0 * rand()) - 1.0d0)
potential_difference = 0.5d0 * (candidate_electric_field_component * candidate_electric_field_component &
                                - electric_field(lattice_site, cartesian_component) &
                                * electric_field(lattice_site, cartesian_component))

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    electric_field(lattice_site, cartesian_component) = candidate_electric_field_component
    no_of_accepted_charge_hops = no_of_accepted_charge_hops + 1
end if

return
end subroutine charge_hop
