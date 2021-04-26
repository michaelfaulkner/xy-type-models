subroutine attempt_external_global_move
use variables
implicit none
integer :: i, lattice_site
double precision :: potential_difference, sign_of_twist
double precision, dimension(no_of_sites) :: candidate_spin_field

potential_difference = 0.0d0
lattice_site = int(dfloat(no_of_sites) * rand()) + 1
sign_of_twist = floor(2.0d0 * rand())
if (sign_of_twist < 0.5d0) then
    sign_of_twist = -1.0d0
end if
! choose a twist in x direction with 0.5 probability
if (floor(2 * rand()) == 0) then
    ! compute and store candidate spin field (with twist applied)
    do i = 1, no_of_sites
        candidate_spin_field(lattice_site) = spin_field(lattice_site) + dfloat(mod(i, integer_lattice_length)) &
                                                                * sign_of_twist * twopi / dfloat(integer_lattice_length)
        lattice_site = lattice_site + mod(lattice_site, no_of_sites) - mod(lattice_site - 1, no_of_sites)
    end do
    ! compute and store candidate emergent-field components and potential difference
    do i = 1, no_of_sites
        potential_difference = potential_difference &
                                        - cos(candidate_spin_field(get_east_neighbour(i)) - candidate_spin_field(i)) &
                                        + cos(spin_field(get_east_neighbour(i)) - spin_field(i))
    end do
    if ((potential_difference < 0.0d0) .or. (rand() < exp(-beta * potential_difference))) then
        do i = 1, no_of_sites
            spin_field(i) = candidate_spin_field(i)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    end if
! else: attempt a twist in y direction
else
    ! compute and store candidate spin field (with twist applied)
    do i = 1, no_of_sites
        candidate_spin_field(lattice_site) = spin_field(lattice_site) &
                + dfloat(mod(int((i - 1) / integer_lattice_length) + 1, integer_lattice_length)) &
                        * sign_of_twist * twopi / dfloat(integer_lattice_length)
        lattice_site = lattice_site + mod(lattice_site, no_of_sites) - mod(lattice_site - 1, no_of_sites)
    end do
    ! compute and store candidate emergent-field components and potential difference
    do i = 1, no_of_sites
        potential_difference = potential_difference &
                                        - cos(candidate_spin_field(get_north_neighbour(i)) - candidate_spin_field(i)) &
                                        + cos(spin_field(get_north_neighbour(i)) - spin_field(i))
    end do
    if ((potential_difference < 0.0d0) .or. (rand() < exp(-beta * potential_difference))) then
        do i = 1, no_of_sites
            spin_field(i) = candidate_spin_field(i)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    end if
end if

return
end subroutine attempt_external_global_move
