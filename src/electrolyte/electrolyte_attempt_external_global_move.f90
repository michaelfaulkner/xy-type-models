subroutine attempt_external_global_move
use variables
implicit none
integer :: i
double precision :: sign_of_topological_sector_change, candidate_electric_field_sum, potential_difference

call calculate_electric_field_sum
potential_difference = 0.0d0
sign_of_topological_sector_change = 2.0d0 * floor(2.0d0 * rand()) - 1.0d0
if (floor(2.0d0 * rand()) == 0) then
    candidate_electric_field_sum = electric_field_sum_x + sign_of_topological_sector_change * elementary_charge &
                                    * integer_lattice_length
    potential_difference = 0.5d0 * (candidate_electric_field_sum * candidate_electric_field_sum &
                                    - electric_field_sum_x * electric_field_sum_x) / dfloat(no_of_sites)
    if ((potential_difference < 0.0d0) .or. (rand() < exp(- beta * potential_difference))) then
        do i = 1, no_of_sites
            electric_field_x(i) = electric_field_x(i) + sign_of_topological_sector_change * elementary_charge &
                                    / integer_lattice_length
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    end if
else
    candidate_electric_field_sum = electric_field_sum_y + sign_of_topological_sector_change * elementary_charge &
                                    * integer_lattice_length
    potential_difference = 0.5d0 * (candidate_electric_field_sum * candidate_electric_field_sum &
                                    - electric_field_sum_y * electric_field_sum_y) / dfloat(no_of_sites)
    if ((potential_difference < 0.0d0) .or. (rand() < exp(- beta * potential_difference))) then
        do i = 1, no_of_sites
           electric_field_y(i) = electric_field_y(i) + sign_of_topological_sector_change * elementary_charge &
                                    / integer_lattice_length
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    end if
end if

return
end subroutine attempt_external_global_move
