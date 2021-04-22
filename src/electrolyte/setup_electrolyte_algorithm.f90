module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: epsilon = 10.0 ** (-6)
integer, parameter :: max_integer_lattice_length = 128
integer, parameter :: max_no_of_sites = max_integer_lattice_length * max_integer_lattice_length
integer :: pos_x(max_no_of_sites), neg_x(max_no_of_sites), pos_y(max_no_of_sites), neg_y(max_no_of_sites)
integer :: rho(max_no_of_sites), array_of_sites(max_no_of_sites)
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps, no_of_observations
integer :: no_of_accepted_field_rotations, no_of_accepted_charge_hops, no_of_accepted_external_global_moves
integer :: ratio_charge_updates, ratio_TSF_updates
double precision :: electric_field_x(max_no_of_sites), electric_field_y(max_no_of_sites)
double precision :: electric_field_sum_x, electric_field_sum_y, elementary_charge
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: volume, width_of_proposal_interval, target_acceptance_rate_of_field_rotations
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
volume = dfloat(no_of_sites)

if (integer_lattice_length > max_integer_lattice_length) then
   write(6,*) 'Linear lattice length exceeds maximum: change the maximum in the common file.'
   stop
end if

return
end subroutine read_in_config_file
