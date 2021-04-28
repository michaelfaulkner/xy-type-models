subroutine attempt_charge_hop(lattice_site, cartesian_component)
use variables
implicit none
integer :: lattice_site, cartesian_component, sign_of_proposed_field_increment
integer :: proposed_charge_value, proposed_neighbouring_charge_value
double precision :: potential_difference, candidate_electric_field_component

sign_of_proposed_field_increment = 2 * int(floor(2.0d0 * rand())) - 1
proposed_charge_value = charge_configuration(lattice_site) - sign_of_proposed_field_increment
if (abs(proposed_charge_value) > 1) then
    return
end if
if (cartesian_component == 1) then
    proposed_neighbouring_charge_value = charge_configuration(get_east_neighbour(lattice_site)) &
                                            + sign_of_proposed_field_increment
else
    proposed_neighbouring_charge_value = charge_configuration(get_north_neighbour(lattice_site)) &
                                            + sign_of_proposed_field_increment
end if
if (abs(proposed_neighbouring_charge_value) > 1) then
    return
end if

candidate_electric_field_component = electric_field(lattice_site, cartesian_component) + elementary_charge &
                                        * sign_of_proposed_field_increment
potential_difference = 0.5d0 * (candidate_electric_field_component * candidate_electric_field_component &
                                - electric_field(lattice_site, cartesian_component) &
                                * electric_field(lattice_site, cartesian_component))

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    electric_field(lattice_site, cartesian_component) = candidate_electric_field_component
    net_charge_displacement(cartesian_component) = net_charge_displacement(cartesian_component) &
                                                    - sign_of_proposed_field_increment
    charge_configuration(lattice_site) = proposed_charge_value
    if (cartesian_component == 1) then
        charge_configuration(get_east_neighbour(lattice_site)) = proposed_neighbouring_charge_value
    else
        charge_configuration(get_north_neighbour(lattice_site)) = proposed_neighbouring_charge_value
    end if
    no_of_accepted_charge_hops = no_of_accepted_charge_hops + 1
end if

return
end subroutine attempt_charge_hop
