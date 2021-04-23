subroutine attempt_external_global_move
use variables
implicit none
integer :: i, lattice_site, initial_lattice_site
double precision :: potential_difference, sign_of_twist

potential_difference = 0.0d0
! choose a twist in x direction with 0.5 probability
if (floor(2 * rand()) == 0) then
    sign_of_twist = floor(2.0d0 * rand())
    if (sign_of_twist < 0.5d0) then
        sign_of_twist = -1.0d0
    end if
    ! pick initial lattice site from which to twist
    initial_lattice_site = int(dfloat(integer_lattice_length) * rand())
    lattice_site = initial_lattice_site
    do i = 1, no_of_sites
        potential_difference = potential_difference - cos(spin_field(get_east_neighbour(lattice_site)) + sign_of_twist &
                                                          * twopi / integer_lattice_length - spin_field(lattice_site)) &
                                                    + cos(spin_field(get_east_neighbour(lattice_site)) &
                                                          - spin_field(lattice_site))
        lattice_site = mod(lattice_site + 1, no_of_sites)
    end do
    if ((potential_difference < 0.0d0) .or. (rand() < exp(-beta * potential_difference))) then
        lattice_site = initial_lattice_site
        do i = 1, no_of_sites
            spin_field(get_east_neighbour(lattice_site)) = mod(spin_field(get_east_neighbour(lattice_site)) &
                                                                + sign_of_twist * twopi / integer_lattice_length, twopi)
            lattice_site = mod(lattice_site + 1, no_of_sites)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    end if
! else: attempt a twist in y direction
else
    sign_of_twist = floor(2.0d0 * rand())
    if (sign_of_twist < 0.5d0) then
        sign_of_twist = -1.0d0
    end if
    ! pick initial lattice site from which to twist
    initial_lattice_site = int(dfloat(integer_lattice_length) * rand())
    lattice_site = initial_lattice_site
    do i = 1, no_of_sites
        potential_difference = potential_difference - cos(spin_field(get_north_neighbour(lattice_site)) + sign_of_twist &
                                                          * twopi / integer_lattice_length - spin_field(lattice_site)) &
                                                    + cos(spin_field(get_north_neighbour(lattice_site)) &
                                                          - spin_field(lattice_site))
        lattice_site = mod(lattice_site + 1, no_of_sites)
    end do
    if ((potential_difference < 0.0d0) .or. (rand() < exp(-beta * potential_difference))) then
        lattice_site = initial_lattice_site
        do i = 1, no_of_sites
            spin_field(get_north_neighbour(lattice_site)) = mod(spin_field(get_north_neighbour(lattice_site)) &
                                                                + sign_of_twist * twopi / integer_lattice_length, twopi)
            lattice_site = mod(lattice_site + 1, no_of_sites)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
    end if
end if

return
end subroutine attempt_external_global_move
