subroutine setup_additional_periodic_boundary
use variables
implicit none
integer :: i

integer_lattice_length = 4
no_of_sites = integer_lattice_length ** 3

! mod(i + integer_lattice_length + 1, integer_lattice_length) rather than mod(i + 1, integer_lattice_length) below as
! mod(i, integer_lattice_length) does not return value in [0, integer_lattice_length) if i < 0
! the up-down-east-west (U-D) notation is defined relative to the following lattice geometry on a 3x3x3 lattice:
! 23 5 14 : D C U (C = current lattice site)

do i = 1, no_of_sites
    get_up_neighbour(i) = i + &
                (mod(int((i - 1) / integer_lattice_length ** 2) + 1, integer_lattice_length) - &
                 mod(int((i - 1) / integer_lattice_length ** 2), integer_lattice_length)) * integer_lattice_length ** 2
    get_down_neighbour(i) = i + &
         (mod(int((i - 1) / integer_lattice_length ** 2) + integer_lattice_length ** 2 - 1, integer_lattice_length) - &
          mod(int((i - 1) / integer_lattice_length ** 2) + integer_lattice_length ** 2, integer_lattice_length)) * &
                integer_lattice_length ** 2
end do

return
end subroutine setup_additional_periodic_boundary
