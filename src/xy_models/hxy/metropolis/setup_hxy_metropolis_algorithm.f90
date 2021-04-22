module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves, randomise_initial_field_configuration, calculate_external_minimising_twist_field
integer, parameter :: max_integer_lattice_length = 128
integer, parameter :: max_no_of_sites = max_integer_lattice_length * max_integer_lattice_length
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: epsilon = 0.00000000001
integer :: pos_x(max_no_of_sites), neg_x(max_no_of_sites), pos_y(max_no_of_sites), neg_y(max_no_of_sites)
integer :: array_of_sites(max_no_of_sites)
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, no_of_accepted_field_rotations, no_of_accepted_external_global_moves, nmax
integer :: no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y
double precision :: theta(max_no_of_sites), emergent_field_x(max_no_of_sites), emergent_field_y(max_no_of_sites)
double precision :: sum_of_squared_electric_field_x, sum_of_squared_electric_field_y
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: volume, width_of_proposal_interval, target_acceptance_rate_of_field_rotations
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
read(10, *) nmax
read(10, *) randomise_initial_field_configuration
read(10, *) use_external_global_moves
read(10, *) calculate_external_minimising_twist_field

if (algorithm_name /= 'hxy-metropolis') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal hxy-metropolis.'
   stop
end if

no_of_sites = integer_lattice_length * integer_lattice_length
volume = dfloat(no_of_sites)
do i = 1, no_of_sites
    array_of_sites(i) = i
end do

if (integer_lattice_length > max_integer_lattice_length) then
    write(6, *) 'Linear lattice length exceeds maximum: change the maximum in the common file.'
    stop
end if

return
end subroutine read_in_config_file
