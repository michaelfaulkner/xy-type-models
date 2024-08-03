module variables
character(100) :: output_directory, algorithm_name, checkpoint_filename
logical :: use_external_global_moves, randomise_initial_field_configuration, measure_magnetisation, measure_helicity
logical :: measure_emergent_field, measure_potential, measure_potential_minimising_twists, measure_external_global_moves
logical :: measure_twist_relaxations, start_from_checkpoint, simulation_complete, print_samples
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: get_up_neighbour, get_down_neighbour, array_of_sites
integer(kind=8) :: no_of_samples ! kind=8 to avoid upper integer bound on long timescales
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: initial_temperature_index, initial_sample_index, external_global_moves(2)
integer :: vacuum_permittivity_sum_cutoff
! kind=8 below to avoid upper integer bound on long timescales
integer(kind=8) :: potential_minimising_twists_squared_sum, potential_minimising_twists_quartic_sum
integer(kind=8) :: twist_relaxations_squared_sum, twist_relaxations_quartic_sum, topological_sector_squared_sum
integer(kind=8) :: topological_sector_quartic_sum
integer :: potential_minimising_twists_sum(2), twist_relaxations_sum(2), topological_sector_sum(2)
real :: time_between_checkpoints, previous_checkpointing_time
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, allocatable, dimension(:) :: spin_field
double precision, allocatable, dimension(:, :) :: emergent_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations, charge_hop_proportion
double precision :: no_of_accepted_field_rotations_per_site, no_of_accepted_external_global_moves
double precision :: raw_magnetic_norm_sum, raw_magnetic_norm_squared_sum, raw_magnetic_norm_quartic_sum
double precision :: raw_inverse_vacuum_perm_sum, raw_inverse_vacuum_perm_squared_sum, raw_macro_josephson_current_sum(2)
double precision :: raw_macro_josephson_current_squared_sum, raw_macro_josephson_current_quartic_sum, potential_sum
double precision :: potential_squared_sum, sum_of_emergent_field_sum(2), sum_of_emergent_field_squared_sum
double precision :: sum_of_emergent_field_quartic_sum
end module variables
