subroutine get_and_print_observation
use variables
implicit none
integer :: i
double precision :: potential

if (measure_electric_field_sum) then
    write(30, 100) - two_pi * dfloat(net_charge_displacement(1)), - two_pi * dfloat(net_charge_displacement(2))
end if

if (measure_potential) then
    potential = 0.0d0
    do i = 1, no_of_sites
        potential = potential + 0.5d0 * (electric_field(i, 1) ** 2 + electric_field(i, 2) ** 2)
    end do
    write(40, 200) potential
end if

if (measure_external_global_moves) then
    write(50, 300) external_global_moves(1), external_global_moves(2)
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)
300 format(I10, ",", I10)

return
end subroutine get_and_print_observation
