subroutine attempt_external_global_move
use variables
implicit none
integer :: i, cartesian_component, sign_of_topological_sector_change, candidate_net_charge_displacement
double precision :: candidate_electric_field_sum, potential_difference

potential_difference = 0.0d0
cartesian_component = 2 * floor(2.0d0 * rand()) - 1
sign_of_topological_sector_change = 2 * floor(2.0d0 * rand()) - 1
candidate_net_charge_displacement = net_charge_displacement(cartesian_component) - sign_of_topological_sector_change &
                                        * integer_lattice_length
potential_difference = 0.5d0 * elementary_charge * elementary_charge &
                        * (candidate_net_charge_displacement * candidate_net_charge_displacement &
                        - net_charge_displacement(cartesian_component) * net_charge_displacement(cartesian_component)) &
                        / dfloat(no_of_sites)

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    do i = 1, no_of_sites
        electric_field(i, cartesian_component) = electric_field(i, cartesian_component) &
                + sign_of_topological_sector_change * elementary_charge / dfloat(integer_lattice_length)
    end do
    net_charge_displacement(cartesian_component) = candidate_net_charge_displacement
    no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
end if

return
end subroutine attempt_external_global_move
