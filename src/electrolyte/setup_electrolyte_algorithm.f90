module variables
character(100) :: output_directory, algorithm_name
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: epsilon = 10.0 ** (-6)
integer, parameter :: max_side = 128
integer, parameter :: max_sites = max_side * max_side
integer :: pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), rho(max_sites), array_of_sites(max_sites)
integer :: side, sites, no_of_temperature_increments, therm_sweeps, measurements, twist
integer :: no_of_accepted_field_rotations, no_of_accepted_charge_hops, no_of_accepted_external_global_moves
integer :: ratio_charge_updates, ratio_TSF_updates
double precision :: Efield_x(max_sites), Efield_y(max_sites), Esum_x, Esum_y, elementaryCharge
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: volume, length, width_of_proposal_interval, target_acceptance_rate_of_field_rotations
end module variables


subroutine read_in_config_file
use variables
implicit none

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) side
read(10, *) therm_sweeps
read(10, *) measurements
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) width_of_proposal_interval
read(10, *) target_acceptance_rate_of_field_rotations
read(10, *) ratio_charge_updates
read(10, *) ratio_TSF_updates
read(10, *) twist
read(10, *) elementaryCharge

if ((algorithm_name /= 'elementary-electrolyte').and.(algorithm_name /= 'multivalued-electrolyte')) then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal either elementary-electrolyte or &
                    multivalued-electrolyte.'
   stop
end if

sites = side * side
length = float(side)
volume = float(sites)

if (side > max_side) then
   write(6,*) 'Linear lattice length exceeds maximum: change the maximum in the common file.'
   stop
end if

return
end subroutine read_in_config_file
