subroutine pre_simulation_processes
use variables
implicit none
character(100) :: config_file, seed_character_string
integer :: seed

! verify that the something has been passed to the exectuable
if (command_argument_count() /= 2) then
    write(6, '(A)') 'Error: pass a configuration file and nine-digit integer-valued seed to executable'
    stop
end if

! read in the seed as a character called seed_character_string
call get_command_argument(1, seed_character_string)
! convert seed_character_string to an integer called seed
read(seed_character_string, *) seed
! setup random-number sequence
call randinit(seed)
write(6, '(A, F16.14)') 'Initial random number = ', rand(seed)

! read in config file
call get_command_argument(2, config_file)
open(unit=10, file=config_file)
call read_config_file

! setup model and lattice
temperature = initial_temperature
if (no_of_temperature_increments == 0) then
    magnitude_of_temperature_increments = 0.0
else
    magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
end if
call setup_periodic_boundaries
call initialise_field_configuration(.true.)

! opens new directory in which to save the sample files
call system('mkdir -p ' // trim(output_directory))

! message informing the start of the Markov process
write(6, '(A, A, A, I4.4, A, I4.4, A)') 'Starting the ', trim(algorithm_name), ' simulation on ', &
                                            integer_lattice_length, 'x', integer_lattice_length, ' lattice sites.'

return
end subroutine pre_simulation_processes