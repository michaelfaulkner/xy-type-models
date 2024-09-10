module variables
character(100) :: output_directory, algorithm_name, checkpoint_filename
logical :: use_external_global_moves, always_hot_start, always_cold_start, measure_magnetisation, measure_helicity
logical :: measure_potential, measure_hot_twist_relaxations, measure_external_global_moves
logical :: measure_twist_relaxations, measure_emergent_field, start_from_checkpoint, simulation_complete, print_samples
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: get_up_neighbour, get_down_neighbour
integer(kind=8) :: no_of_samples ! kind=8 to avoid upper integer bound on long timescales
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: initial_temperature_index, initial_sample_index, external_global_moves(2)
! kind=8 below to avoid upper integer bound on long timescales
integer(kind=8) :: hot_twist_relaxations_squared_sum, hot_twist_relaxations_quartic_sum, twist_relaxations_squared_sum
integer(kind=8) :: twist_relaxations_quartic_sum, topological_sector_squared_sum, topological_sector_quartic_sum
integer :: hot_twist_relaxations_sum(2), twist_relaxations_sum(2), topological_sector_sum(2)
real :: time_between_checkpoints, previous_checkpointing_time
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, allocatable, dimension(:) :: spin_field
double precision, allocatable, dimension(:, :) :: emergent_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: spin_space_distance_between_samples, no_of_events_per_unit_spin_space_distance
double precision :: no_of_accepted_external_global_moves
double precision :: magnetic_norm_sum, magnetic_norm_squared_sum, magnetic_norm_quartic_sum
double precision :: inverse_vacuum_perm_sum, inverse_vacuum_perm_squared_sum, macro_josephson_current_sum(2)
double precision :: macro_josephson_current_squared_sum, macro_josephson_current_quartic_sum, potential_sum
double precision :: potential_squared_sum, potential_quartic_sum, emergent_field_zero_mode_sum(2)
double precision :: emergent_field_zero_mode_squared_sum, emergent_field_zero_mode_quartic_sum
end module variables
