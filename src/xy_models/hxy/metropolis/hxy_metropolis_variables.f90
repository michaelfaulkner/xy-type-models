module variables
character(100) :: output_directory, algorithm_name, checkpoint_filename
logical :: use_external_global_moves, randomise_initial_field_configuration, measure_magnetisation, measure_helicity
logical :: measure_potential, measure_potential_minimising_twists, measure_external_global_moves, start_from_checkpoint
logical :: simulation_complete
integer, allocatable, dimension(:) :: get_north_neighbour, get_south_neighbour, get_east_neighbour, get_west_neighbour
integer, allocatable, dimension(:) :: get_up_neighbour, get_down_neighbour, array_of_sites
integer(kind=8) :: no_of_observations ! kind=8 to avoid upper integer bound on long timescales
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: initial_temperature_index, initial_observation_index, external_global_moves(2)
integer :: vacuum_permittivity_sum_cutoff
real :: time_between_checkpoints, previous_checkpointing_time
double precision, parameter :: two_pi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, allocatable, dimension(:) :: spin_field
double precision, allocatable, dimension(:, :) :: emergent_field
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations, charge_hop_proportion
double precision :: no_of_accepted_field_rotations_per_site, no_of_accepted_external_global_moves
end module variables
