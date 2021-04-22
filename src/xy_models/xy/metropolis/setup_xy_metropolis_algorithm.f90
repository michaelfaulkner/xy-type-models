module variables
integer, allocatable, dimension(:) :: pos_x, neg_x, pos_y, neg_y, array_of_sites
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, no_of_accepted_field_rotations, no_of_accepted_external_global_moves
double precision, parameter :: twopi = 6.28318530717959d0
double precision, allocatable, dimension(:) :: theta
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves, randomise_initial_field_configuration
end module variables


subroutine read_in_config_file
use variables
implicit none
integer :: i

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) integer_lattice_length
read(10, *) no_of_equilibration_sweeps
read(10, *) no_of_observations
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) width_of_proposal_interval
read(10, *) target_acceptance_rate_of_field_rotations
read(10, *) randomise_initial_field_configuration
read(10, *) use_external_global_moves

if (algorithm_name /= 'xy-metropolis') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal xy-metropolis.'
   stop
end if

no_of_sites = integer_lattice_length * integer_lattice_length
allocate(theta(no_of_sites), pos_x(no_of_sites), pos_y(no_of_sites), neg_x(no_of_sites), neg_y(no_of_sites), &
            array_of_sites(no_of_sites))
do i = 1, no_of_sites
    array_of_sites(i) = i
end do

return
end subroutine read_in_config_file
