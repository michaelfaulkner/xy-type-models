subroutine attempt_single_external_global_move(cartesian_component)
use variables
implicit none
integer :: cartesian_component, i, lattice_site, sign_of_twist
double precision :: potential_difference, get_spin_difference
double precision, dimension(no_of_sites) :: candidate_spin_field, candidate_emergent_field_components

potential_difference = 0.0d0
sign_of_twist = 2 * floor(2.0d0 * rand()) - 1
if (cartesian_component == 1) then
    ! attempt a global twist along x direction
    lattice_site = int(dfloat(no_of_sites) * rand()) + 1
    ! compute and store candidate spin field (with twist applied)
    do i = 1, no_of_sites
        candidate_spin_field(lattice_site) = spin_field(lattice_site) + dfloat(mod(i, integer_lattice_length)) &
                                                    * dfloat(sign_of_twist) * two_pi / dfloat(integer_lattice_length)
        lattice_site = lattice_site + mod(lattice_site, no_of_sites) - mod(lattice_site - 1, no_of_sites)
    end do
    ! compute and store candidate emergent-field components and potential difference
    do i = 1, no_of_sites
        candidate_emergent_field_components(i) = get_spin_difference(candidate_spin_field(get_west_neighbour(i)), &
                                                                     candidate_spin_field(i))
        potential_difference = potential_difference + 0.5d0 * &
                                            (candidate_emergent_field_components(i) ** 2 - emergent_field(i, 2) ** 2)
    end do
else if (cartesian_component == 2) then
    ! attempt a global twist along y direction
    lattice_site = floor(dfloat(no_of_sites) * rand() / dfloat(integer_lattice_length)) * integer_lattice_length + 1
    ! compute and store candidate spin field (with twist applied)
    do i = 1, no_of_sites
        candidate_spin_field(lattice_site) = spin_field(lattice_site) &
                                    + dfloat(mod(int((i - 1) / integer_lattice_length) + 1, integer_lattice_length)) &
                                            * dfloat(sign_of_twist) * two_pi / dfloat(integer_lattice_length)
        lattice_site = lattice_site + mod(lattice_site, no_of_sites) - mod(lattice_site - 1, no_of_sites)
    end do
    ! compute and store candidate emergent-field components and potential difference
    do i = 1, no_of_sites
        candidate_emergent_field_components(i) = get_spin_difference(candidate_spin_field(i), &
                                                                     candidate_spin_field(get_south_neighbour(i)))
        potential_difference = potential_difference + 0.5d0 * &
                                            (candidate_emergent_field_components(i) ** 2 - emergent_field(i, 1) ** 2)
    end do
end if

if ((potential_difference < 0.0d0).or.(rand() < exp(-beta * potential_difference))) then
    do i = 1, no_of_sites
        spin_field(i) = candidate_spin_field(i)
        emergent_field(i, 3 - cartesian_component) = candidate_emergent_field_components(i)
    end do
    ! we count accepted Metropolis moves in double precision (float) to avoid upper integer bound on long timescales
    no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1.0d0
    external_global_moves(cartesian_component) = sign_of_twist
else
    external_global_moves(cartesian_component) = 0
end if

return
end subroutine attempt_single_external_global_move
