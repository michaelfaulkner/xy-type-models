subroutine metropolis_sweep
use variables
implicit none
integer :: n, i
double precision :: candidate_spin_value, potential_difference, get_observation_of_normal_distribution

call randomise_array_of_sites
do n = 1, no_of_sites
    i = array_of_sites(n)
    candidate_spin_value = mod(spin_field(i) + &
                                get_observation_of_normal_distribution(0.0d0, width_of_proposal_interval), two_pi)

    potential_difference = - cos(spin_field(get_north_neighbour(i)) - candidate_spin_value) &
                            - cos(candidate_spin_value - spin_field(get_south_neighbour(i))) &
                            - cos(spin_field(get_east_neighbour(i)) - candidate_spin_value) &
                            - cos(candidate_spin_value - spin_field(get_west_neighbour(i))) &
                            - cos(spin_field(get_up_neighbour(i)) - candidate_spin_value) &
                            - cos(candidate_spin_value - spin_field(get_down_neighbour(i))) &
                            + cos(spin_field(get_north_neighbour(i)) - spin_field(i)) &
                            + cos(spin_field(i) - spin_field(get_south_neighbour(i))) &
                            + cos(spin_field(get_east_neighbour(i)) - spin_field(i)) &
                            + cos(spin_field(i) - spin_field(get_west_neighbour(i))) &
                            + cos(spin_field(get_up_neighbour(i)) - spin_field(i)) &
                            + cos(spin_field(i) - spin_field(get_down_neighbour(i)))

    if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
        spin_field(i) = candidate_spin_value
        no_of_accepted_field_rotations = no_of_accepted_field_rotations + 1
    end if
end do

return
end subroutine metropolis_sweep
