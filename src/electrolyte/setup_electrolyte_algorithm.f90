module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves
integer, allocatable, dimension(:) :: pos_x, neg_x, pos_y, neg_y, array_of_sites, rho
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, no_of_accepted_field_rotations, no_of_accepted_external_global_moves
integer :: no_of_accepted_charge_hops, ratio_charge_updates, ratio_TSF_updates
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, allocatable, dimension(:) :: electric_field_x, electric_field_y
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: width_of_proposal_interval, target_acceptance_rate_of_field_rotations
double precision :: electric_field_sum_x, electric_field_sum_y, elementary_charge
end module variables


subroutine read_in_config_file
use variables
implicit none

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
read(10, *) ratio_charge_updates
read(10, *) ratio_TSF_updates
read(10, *) use_external_global_moves
read(10, *) elementary_charge

if ((algorithm_name /= 'elementary-electrolyte').and.(algorithm_name /= 'multivalued-electrolyte')) then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal either elementary-electrolyte or &
                    multivalued-electrolyte.'
   stop
end if

no_of_sites = integer_lattice_length * integer_lattice_length
allocate(electric_field_x(no_of_sites), electric_field_y(no_of_sites))
allocate(pos_x(no_of_sites), pos_y(no_of_sites), neg_x(no_of_sites), neg_y(no_of_sites), array_of_sites(no_of_sites))
allocate(rho(no_of_sites))

return
end subroutine read_in_config_file
