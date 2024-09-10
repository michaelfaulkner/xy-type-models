module variables
character(100) :: output_directory, algorithm_name, checkpoint_filename
logical :: always_cold_start, use_external_global_moves, measure_electric_field_sum, measure_potential
logical :: measure_external_global_moves, start_from_checkpoint, simulation_complete, print_samples
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: get_up_neighbour, get_down_neighbour, array_of_sites, charge_configuration
integer(kind=8) :: no_of_samples ! kind=8 to avoid upper integer bound on long timescales
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: initial_temperature_index, initial_sample_index, net_charge_displacement(2), external_global_moves(2)
integer :: topological_sector_sum(2)
! kind=8 below to avoid upper integer bound on long timescales
integer(kind=8) :: topological_sector_squared_sum, topological_sector_quartic_sum
real :: time_between_checkpoints, previous_checkpointing_time
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, allocatable, dimension(:, :) :: electric_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations
double precision :: charge_hop_proportion, charge_hop_proportion_over_two, no_of_accepted_charge_hops_per_site
double precision :: no_of_accepted_field_rotations_per_site, no_of_accepted_external_global_moves
double precision :: potential_sum, potential_squared_sum, potential_quartic_sum, raw_electric_field_zero_mode_sum(2)
double precision :: raw_electric_field_zero_mode_squared_sum, raw_electric_field_zero_mode_quartic_sum
end module variables
