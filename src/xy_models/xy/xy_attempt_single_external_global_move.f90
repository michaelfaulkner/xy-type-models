subroutine attempt_single_external_global_move(cartesian_component)
use variables
implicit none
integer :: cartesian_component, i, lattice_site, sign_of_twist
double precision :: potential_difference
double precision, dimension(no_of_sites) :: candidate_spin_field

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
        potential_difference = potential_difference &
                                        - cos(candidate_spin_field(get_west_neighbour(i)) - candidate_spin_field(i)) &
                                        + cos(spin_field(get_west_neighbour(i)) - spin_field(i))
    end do
    if ((potential_difference < 0.0d0).or.(rand() < exp(-beta * potential_difference))) then
        do i = 1, no_of_sites
            spin_field(i) = candidate_spin_field(i)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
        external_global_move(1) = sign_of_twist
    else
        external_global_move(1) = 0
    end if
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
        potential_difference = potential_difference &
                                        - cos(candidate_spin_field(i) - candidate_spin_field(get_south_neighbour(i))) &
                                        + cos(spin_field(i) - spin_field(get_south_neighbour(i)))
    end do
    if ((potential_difference < 0.0d0).or.(rand() < exp(-beta * potential_difference))) then
        do i = 1, no_of_sites
            spin_field(i) = candidate_spin_field(i)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
        external_global_move(2) = sign_of_twist
    else
        external_global_move(2) = 0
    end if
end if

return
end subroutine attempt_single_external_global_move
