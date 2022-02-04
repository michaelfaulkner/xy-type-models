subroutine attempt_single_external_global_move(cartesian_component)
use variables
implicit none
integer :: cartesian_component, i, sign_of_topological_sector_change, candidate_net_charge_displacement
double precision :: candidate_electric_field_sum, potential_difference

potential_difference = 0.0d0
sign_of_topological_sector_change = 2 * floor(2.0d0 * rand()) - 1
candidate_net_charge_displacement = net_charge_displacement(cartesian_component) - sign_of_topological_sector_change &
                                        * integer_lattice_length
potential_difference = 0.5d0 * two_pi * two_pi * (dfloat(candidate_net_charge_displacement) ** 2 &
                        - dfloat(net_charge_displacement(cartesian_component)) ** 2) / dfloat(no_of_sites)

if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
    do i = 1, no_of_sites
        electric_field(i, cartesian_component) = electric_field(i, cartesian_component) &
                + dfloat(sign_of_topological_sector_change) * two_pi / dfloat(integer_lattice_length)
    end do
    net_charge_displacement(cartesian_component) = candidate_net_charge_displacement
    no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    external_global_move(cartesian_component) = sign_of_topological_sector_change
else
    external_global_move(cartesian_component) = 0
end if

return
end subroutine attempt_single_external_global_move
