subroutine setup_periodic_boundaries
use variables
implicit none
integer :: i

! mod(i + integer_lattice_length + 1, integer_lattice_length) rather than mod(i + 1, integer_lattice_length) below as
! mod(i, integer_lattice_length) does not return value in [0, integer_lattice_length) if i < 0
do i = 1, no_of_sites
    pos_x(i) = i + mod(i, integer_lattice_length) - mod(i - 1, integer_lattice_length)
    neg_x(i) = i + mod(i - 2 + integer_lattice_length, integer_lattice_length) - &
                    mod(i - 1 + integer_lattice_length, integer_lattice_length)
    pos_y(i) = i + (mod(int((i - 1) / integer_lattice_length) + 1, integer_lattice_length) - &
                    mod(int((i - 1) / integer_lattice_length), integer_lattice_length)) * integer_lattice_length
    neg_y(i) = i + (mod(int((i - 1) / integer_lattice_length) + integer_lattice_length - 1, integer_lattice_length) - &
                    mod(int((i - 1) / integer_lattice_length) + integer_lattice_length, integer_lattice_length)) * &
                    integer_lattice_length
end do

return
end subroutine setup_periodic_boundaries
