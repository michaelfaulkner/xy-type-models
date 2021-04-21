subroutine pre_simulation_processes
use variables
implicit none
character(100) :: config_file
integer :: seed_size
integer, allocatable :: seed(:)

! verify that the something has been parsed to the exectuable
if (command_argument_count() /= 1) then
    write(6, *) 'Error: parse configuration file to executable'
    stop
end if

! setup random-number sequence
seed_size = 123456
call random_seed(size=seed_size)
allocate(seed(seed_size))
call random_seed(get=seed)
call randinit(abs(seed(1)))
write(6, '(A, F16.14)') 'Initial random number = ', rand(abs(seed(1)))

! read in config file
call get_command_argument(1, config_file)
open (unit=1, file=config_file)
call read_in_config_file

! setup model and lattice
temperature = initial_temperature
if (no_of_temperature_increments == 0) then
    magnitude_of_temperature_increments = 0.0
else
    magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
end if
call setup_periodic_boundaries
call initialise_field_configuration

return
end subroutine pre_simulation_processes