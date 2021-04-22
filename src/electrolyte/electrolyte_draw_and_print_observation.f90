subroutine draw_and_print_observation
use variables
implicit none
integer i
double precision potential

potential = 0.0d0
do i = 1, no_of_sites
    potential = potential + 0.5d0 * (electric_field_x(i) * electric_field_x(i) + &
                                            electric_field_y(i) * electric_field_y(i))
end do

call calculate_electric_field_sum

write(20, 100) potential, electric_field_sum_x, electric_field_sum_y

100 format(ES24.14, ", ", ES24.14, ", ", ES24.14)

return
end subroutine draw_and_print_observation

subroutine calculate_electric_field_sum
use variables
implicit none
integer i

electric_field_sum_x = 0.0d0
electric_field_sum_y = 0.0d0
do i = 1, no_of_sites
    electric_field_sum_x = electric_field_sum_x + electric_field_x(i)
    electric_field_sum_y = electric_field_sum_y + electric_field_y(i)
end do

return
end subroutine calculate_electric_field_sum
