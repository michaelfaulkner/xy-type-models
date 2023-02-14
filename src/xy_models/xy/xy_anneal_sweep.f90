subroutine xy_anneal_sweep
use variables
implicit none
integer :: i
double precision :: candidate_spin_value, potential_difference, get_potential_difference, get_spin_difference
double precision :: spin_perturbation, initial_spin_difference_1, initial_spin_difference_2, initial_spin_difference_3
double precision :: initial_spin_difference_4, perturbed_spin_difference_1, perturbed_spin_difference_2
double precision :: perturbed_spin_difference_3, perturbed_spin_difference_4

do i = 1, no_of_sites
    spin_perturbation = 0.01d0 * (rand() - 0.5d0)

    initial_spin_difference_1 = get_spin_difference(spin_field(i), spin_field(get_south_neighbour(i)))
    initial_spin_difference_2 = get_spin_difference(spin_field(get_west_neighbour(i)), spin_field(i))
    initial_spin_difference_3 = get_spin_difference(spin_field(get_north_neighbour(i)), spin_field(i))
    initial_spin_difference_4 = get_spin_difference(spin_field(i), spin_field(get_east_neighbour(i)))

    perturbed_spin_difference_1 = initial_spin_difference_1 + spin_perturbation
    perturbed_spin_difference_2 = initial_spin_difference_2 - spin_perturbation
    perturbed_spin_difference_3 = initial_spin_difference_3 - spin_perturbation
    perturbed_spin_difference_4 = initial_spin_difference_4 + spin_perturbation

    if ((abs(perturbed_spin_difference_1) < pi).and.(abs(perturbed_spin_difference_2) < pi).and.&
            (abs(perturbed_spin_difference_3) < pi).and.(abs(perturbed_spin_difference_4) < pi)) then
        ! we may continue as candidate move does not alter the topological-defect configuration
        candidate_spin_value = mod(spin_field(i) + spin_perturbation, two_pi)
        potential_difference = get_potential_difference(candidate_spin_value, i)
        if (potential_difference < 0.0d0) then
            spin_field(i) = candidate_spin_value
        end if
    end if
end do

return
end subroutine xy_anneal_sweep
