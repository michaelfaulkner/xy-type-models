subroutine metropolis_sweep
use variables
implicit none
integer :: n, i
double precision :: candidate_theta, potential_difference, candidate_emergent_field_1, candidate_emergent_field_2
double precision :: candidate_emergent_field_3, candidate_emergent_field_4

call randomise_array_of_sites
do n = 1, sites
    i = array_of_sites(n)
    candidate_theta = mod(theta(i) + width_of_proposal_interval * (rand() - 0.5d0), twopi)
    candidate_emergent_field_1 = modulo(candidate_theta - theta(neg_y(i)) + pi, twopi) - pi
    candidate_emergent_field_2 = modulo(- candidate_theta + theta(neg_x(i)) + pi, twopi) - pi
    candidate_emergent_field_3 = modulo(theta(pos_y(i)) - candidate_theta + pi, twopi) - pi
    candidate_emergent_field_4 = modulo(- theta(pos_x(i)) + candidate_theta + pi, twopi) - pi

    potential_difference = 0.5d0 * (candidate_emergent_field_1 * candidate_emergent_field_1 &
                                    + candidate_emergent_field_2 * candidate_emergent_field_2 &
                                    + candidate_emergent_field_3 * candidate_emergent_field_3 &
                                    + candidate_emergent_field_4 * candidate_emergent_field_4 &
                                    - emergent_field_x(i) * emergent_field_x(i) &
                                    - emergent_field_y(i) * emergent_field_y(i) &
                                    - emergent_field_x(pos_y(i)) * emergent_field_x(pos_y(i)) &
                                    - emergent_field_y(pos_x(i)) * emergent_field_y(pos_x(i)))

    if ((potential_difference < 0.0d0) .or. (rand() < exp(- beta * potential_difference))) then
        theta(i) = candidate_theta
        emergent_field_x(i) = candidate_emergent_field_1
        emergent_field_y(i) = candidate_emergent_field_2
        emergent_field_x(pos_y(i)) = candidate_emergent_field_3
        emergent_field_y(pos_x(i)) = candidate_emergent_field_4
        no_of_accepted_field_rotations = no_of_accepted_field_rotations + 1
    end if
end do

return
end subroutine metropolis_sweep
