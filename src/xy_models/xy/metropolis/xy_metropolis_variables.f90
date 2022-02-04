module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves, randomise_initial_field_configuration
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: array_of_sites
integer(kind=8) :: no_of_accepted_field_rotations
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, no_of_accepted_external_global_moves, external_global_move(2)
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, allocatable, dimension(:) :: spin_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations, charge_hop_proportion
end module variables
