module variables
character(100) :: output_directory, algorithm_name
logical :: randomise_initial_field_configuration
logical, parameter :: use_external_global_moves = .false.
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: get_up_neighbour, get_down_neighbour, array_of_sites
integer(kind=8) :: no_of_accepted_field_rotations
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, allocatable, dimension(:) :: spin_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations, charge_hop_proportion
end module variables
