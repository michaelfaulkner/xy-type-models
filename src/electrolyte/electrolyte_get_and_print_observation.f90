subroutine get_and_print_observation
use variables
implicit none
integer :: i
double precision :: potential

potential = 0.0d0
do i = 1, no_of_sites
    potential = potential + 0.5d0 * (electric_field(i, 1) * electric_field(i, 1) &
                                    + electric_field(i, 2) * electric_field(i, 2))
end do

write(20, 100) potential, - elementary_charge * dfloat(net_charge_displacement(1)), &
                          - elementary_charge * dfloat(net_charge_displacement(2))

100 format(ES24.14, ", ", ES24.14, ", ", ES24.14)

return
end subroutine get_and_print_observation
