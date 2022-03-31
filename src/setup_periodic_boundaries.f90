subroutine setup_periodic_boundaries
use variables
implicit none
integer :: i

! mod(i + integer_lattice_length + 1, integer_lattice_length) rather than mod(i + 1, integer_lattice_length) below as
! mod(i, integer_lattice_length) does not return value in [0, integer_lattice_length) if i < 0
! the north-south-east-west (N-S-E-W) notation is defined relative to the following lattice geometry:
! 7 8 9    N
! 4 5 6  W   E
! 1 2 3    S
! the up-down-east-west (U-D) notation is defined relative to the following lattice geometry on a 3x3x3 lattice:
! 23 5 14 : D C U (C = current lattice site)

do i = 1, no_of_sites
    get_east_neighbour(i) = i + mod(i, integer_lattice_length) - mod(i - 1, integer_lattice_length)
    get_west_neighbour(i) = i + mod(i - 2 + integer_lattice_length, integer_lattice_length) - &
                                mod(i - 1 + integer_lattice_length, integer_lattice_length)
    get_north_neighbour(i) = i + &
            (mod(int((i - 1) / integer_lattice_length) + 1, integer_lattice_length) - &
                mod(int((i - 1) / integer_lattice_length), integer_lattice_length)) * integer_lattice_length
    get_south_neighbour(i) = i + &
            (mod(int((i - 1) / integer_lattice_length) + integer_lattice_length - 1, integer_lattice_length) - &
                mod(int((i - 1) / integer_lattice_length) + integer_lattice_length, integer_lattice_length)) * &
                    integer_lattice_length
    get_up_neighbour(i) = i + &
                (mod(int((i - 1) / integer_lattice_length ** 2) + 1, integer_lattice_length) - &
                 mod(int((i - 1) / integer_lattice_length ** 2), integer_lattice_length)) * integer_lattice_length ** 2
    get_down_neighbour(i) = i + &
         (mod(int((i - 1) / integer_lattice_length ** 2) + integer_lattice_length ** 2 - 1, integer_lattice_length) - &
          mod(int((i - 1) / integer_lattice_length ** 2) + integer_lattice_length ** 2, integer_lattice_length)) * &
                integer_lattice_length ** 2
end do

return
end subroutine setup_periodic_boundaries
