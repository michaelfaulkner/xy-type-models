subroutine draw_and_print_observation
use variables
implicit none
integer i
double precision potential

potential = 0.0d0
do i = 1, sites
    potential = potential + 0.5 * (Efield_x(i) * Efield_x(i) + Efield_y(i) * Efield_y(i))
end do

call calculate_electric_field_sum

write(20, 100) potential, Esum_x, Esum_y

100 format(ES24.14, ", ", ES24.14, ", ", ES24.14)

return
end subroutine draw_and_print_observation

subroutine calculate_electric_field_sum
use variables
implicit none
integer i

Esum_x = 0.0d0
Esum_y = 0.0d0
do i = 1, sites
    Esum_x = Esum_x + Efield_x(i)
    Esum_y = Esum_y + Efield_y(i)
end do

return
end subroutine calculate_electric_field_sum
