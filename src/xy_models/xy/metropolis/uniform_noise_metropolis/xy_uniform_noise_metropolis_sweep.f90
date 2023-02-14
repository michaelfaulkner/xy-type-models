subroutine metropolis_sweep
use variables
implicit none
integer :: n, i
double precision :: candidate_spin_value, potential_difference, get_potential_difference

call randomise_array_of_sites
do n = 1, no_of_sites
    i = array_of_sites(n)
    candidate_spin_value = mod(spin_field(i) + width_of_proposal_interval * (rand() - 0.5d0), two_pi)
    potential_difference = get_potential_difference(candidate_spin_value, i)
    if ((potential_difference < 0.0d0).or.(rand() < exp(- beta * potential_difference))) then
        spin_field(i) = candidate_spin_value
        ! we count accepted Metropolis moves in double precision (float) to avoid upper integer bound on long timescales
        no_of_accepted_field_rotations_per_site = no_of_accepted_field_rotations_per_site + 1.0d0 / dfloat(no_of_sites)
    end if
end do

return
end subroutine metropolis_sweep
