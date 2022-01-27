subroutine attempt_charge_hop(lattice_site, cartesian_component)
use variables
implicit none
integer :: lattice_site, cartesian_component, sign_of_proposed_integer_field_increment
double precision :: potential_difference, candidate_electric_field_component

sign_of_proposed_integer_field_increment = 2 * floor(2.0d0 * rand()) - 1
candidate_electric_field_component = electric_field(lattice_site, cartesian_component) + two_pi &
                                        * dfloat(sign_of_proposed_integer_field_increment)
potential_difference = 0.5d0 * (candidate_electric_field_component ** 2 &
                                - electric_field(lattice_site, cartesian_component) ** 2)

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    electric_field(lattice_site, cartesian_component) = candidate_electric_field_component
    net_charge_displacement(cartesian_component) = net_charge_displacement(cartesian_component) &
                                                    - sign_of_proposed_integer_field_increment
    no_of_accepted_charge_hops = no_of_accepted_charge_hops + 1
end if

return
end subroutine attempt_charge_hop
