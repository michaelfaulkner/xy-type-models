function get_potential_difference(candidate_spin_value, site_index)
use variables
implicit none
integer :: site_index
double precision :: get_potential_difference, candidate_spin_value

get_potential_difference = - cos(spin_field(get_east_neighbour(site_index)) - candidate_spin_value) &
                            - cos(spin_field(get_north_neighbour(site_index)) - candidate_spin_value) &
                            - cos(candidate_spin_value - spin_field(get_west_neighbour(site_index))) &
                            - cos(candidate_spin_value - spin_field(get_south_neighbour(site_index))) &
                            + cos(spin_field(get_east_neighbour(site_index)) - spin_field(site_index)) &
                            + cos(spin_field(get_north_neighbour(site_index)) - spin_field(site_index)) &
                            + cos(spin_field(site_index) - spin_field(get_west_neighbour(site_index))) &
                            + cos(spin_field(site_index) - spin_field(get_south_neighbour(site_index)))

end function get_potential_difference
