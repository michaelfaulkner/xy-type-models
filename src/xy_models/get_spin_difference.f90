function get_spin_difference(spin_value_one, spin_value_two)
use variables
implicit none
double precision :: get_spin_difference, spin_value_one, spin_value_two
get_spin_difference = modulo(spin_value_one - spin_value_two + pi, twopi) - pi
end function get_spin_difference
