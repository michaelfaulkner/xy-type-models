subroutine attempt_charge_hop(lattice_site, cartesian_component)
use variables
implicit none
integer :: lattice_site, cartesian_component, sign_of_proposed_integer_field_increment
integer :: proposed_integer_charge_value, proposed_neighbouring_integer_charge_value
double precision :: potential_difference, candidate_electric_field_component

sign_of_proposed_integer_field_increment = 2 * floor(2.0d0 * rand()) - 1
proposed_integer_charge_value = charge_configuration(lattice_site) + sign_of_proposed_integer_field_increment
if (abs(proposed_integer_charge_value) > 1) then
    return
end if
if (cartesian_component == 1) then
    proposed_neighbouring_integer_charge_value = charge_configuration(get_east_neighbour(lattice_site)) &
                                                    - sign_of_proposed_integer_field_increment
else
    proposed_neighbouring_integer_charge_value = charge_configuration(get_north_neighbour(lattice_site)) &
                                                    - sign_of_proposed_integer_field_increment
end if
if (abs(proposed_neighbouring_integer_charge_value) > 1) then
    return
end if

candidate_electric_field_component = electric_field(lattice_site, cartesian_component) + two_pi &
                                        * dfloat(sign_of_proposed_integer_field_increment)
potential_difference = 0.5d0 * (candidate_electric_field_component ** 2 &
                                - electric_field(lattice_site, cartesian_component) ** 2)

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    electric_field(lattice_site, cartesian_component) = candidate_electric_field_component
    net_charge_displacement(cartesian_component) = net_charge_displacement(cartesian_component) &
                                                    - sign_of_proposed_integer_field_increment
    charge_configuration(lattice_site) = proposed_integer_charge_value
    if (cartesian_component == 1) then
        charge_configuration(get_east_neighbour(lattice_site)) = proposed_neighbouring_integer_charge_value
    else
        charge_configuration(get_north_neighbour(lattice_site)) = proposed_neighbouring_integer_charge_value
    end if
    ! we count accepted Metropolis moves in double precision (float) to avoid upper integer bound on long timescales
    no_of_accepted_charge_hops_per_site = no_of_accepted_charge_hops_per_site + 1.0d0 / dfloat(no_of_sites)
end if

return
end subroutine attempt_charge_hop
