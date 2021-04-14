subroutine metropolis_sweep
use variables
implicit none
integer :: n, i
double precision :: candidate_theta, potential_difference

call randomise_array_of_sites
do n = 1, sites
    i = array_of_sites(n)
    candidate_theta = mod(theta(i) + width_of_proposal_interval * (rand() - 0.5d0), twopi)

    potential_difference = - cos(theta(pos_x(i)) - candidate_theta) - cos(theta(pos_y(i)) - candidate_theta) &
                                - cos(candidate_theta - theta(neg_x(i))) - cos(candidate_theta - theta(neg_y(i))) &
                                + cos(theta(pos_x(i)) - theta(i)) + cos(theta(pos_y(i)) - theta(i)) &
                                + cos(theta(i) - theta(neg_x(i))) + cos(theta(i) - theta(neg_y(i)))

    if ((potential_difference < 0.0d0) .or. (rand() < exp(- beta * potential_difference))) then
        theta(i) = candidate_theta
        no_of_accepted_field_rotations = no_of_accepted_field_rotations + 1
    end if
end do

return
end subroutine metropolis_sweep
