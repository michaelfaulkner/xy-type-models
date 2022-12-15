module variables
character(100) :: output_directory, algorithm_name, checkpoint_filename
logical :: use_external_global_moves, measure_electric_field_sum, measure_potential, measure_external_global_moves
logical :: start_from_checkpoint
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: get_up_neighbour, get_down_neighbour, array_of_sites, charge_configuration
integer(kind=8) :: no_of_accepted_charge_hops, no_of_accepted_field_rotations, no_of_accepted_external_global_moves
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, initial_temperature_index, initial_observation_index, net_charge_displacement(2)
integer :: external_global_moves(2)
real :: time_between_checkpoints, previous_checkpointing_time
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, allocatable, dimension(:, :) :: electric_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations
double precision :: charge_hop_proportion, charge_hop_proportion_over_two
end module variables
