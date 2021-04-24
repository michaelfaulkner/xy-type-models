module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves, randomise_initial_field_configuration, calculate_external_minimising_twist_field
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, no_of_events, no_of_accepted_external_global_moves, vacuum_permittivity_sum_cutoff
integer :: no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: pi_squared = 9.86960440108936d0
double precision, parameter :: pi_squared_over_two = 4.93480220054468d0
double precision, allocatable, dimension(:) :: spin_field, emergent_field_x, emergent_field_y
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: spin_space_distance_between_observations
double precision :: sum_of_squared_electric_field_x, sum_of_squared_electric_field_y
end module variables