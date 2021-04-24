subroutine metropolis_sweep
use variables
implicit none
integer :: n, i
double precision :: get_spin_difference, candidate_spin_value, potential_difference, candidate_emergent_field(4)

call randomise_array_of_sites
do n = 1, no_of_sites
    i = array_of_sites(n)
    candidate_spin_value = mod(spin_field(i) + width_of_proposal_interval * (rand() - 0.5d0), twopi)
    candidate_emergent_field(1) = get_spin_difference(candidate_spin_value, spin_field(get_south_neighbour(i)))
    candidate_emergent_field(2) = get_spin_difference(spin_field(get_west_neighbour(i)), candidate_spin_value)
    candidate_emergent_field(3) = get_spin_difference(spin_field(get_north_neighbour(i)), candidate_spin_value)
    candidate_emergent_field(4) = get_spin_difference(candidate_spin_value, spin_field(get_east_neighbour(i)))

    potential_difference = 0.5d0 * (&
                              candidate_emergent_field(1) * candidate_emergent_field(1) &
                            + candidate_emergent_field(2) * candidate_emergent_field(2) &
                            + candidate_emergent_field(3) * candidate_emergent_field(3) &
                            + candidate_emergent_field(4) * candidate_emergent_field(4) &
                            - emergent_field(i, 1) * emergent_field(i, 1) &
                            - emergent_field(i, 2) * emergent_field(i, 2) &
                            - emergent_field(get_north_neighbour(i), 1) * emergent_field(get_north_neighbour(i), 1) &
                            - emergent_field(get_east_neighbour(i), 2) * emergent_field(get_east_neighbour(i), 2))

    if ((potential_difference < 0.0d0) .or. (rand() < exp(- beta * potential_difference))) then
        spin_field(i) = candidate_spin_value
        emergent_field(i, 1) = candidate_emergent_field(1)
        emergent_field(i, 2) = candidate_emergent_field(2)
        emergent_field(get_north_neighbour(i), 1) = candidate_emergent_field(3)
        emergent_field(get_east_neighbour(i), 2) = candidate_emergent_field(4)
        no_of_accepted_field_rotations = no_of_accepted_field_rotations + 1
    end if
end do

return
end subroutine metropolis_sweep
