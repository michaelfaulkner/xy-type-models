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

! search for checkpoint
time_between_checkpoints = 1.0
previous_checkpointing_time = 0.0
write(checkpoint_filename, '(A, "/checkpoint.csv")') trim(output_directory)
inquire(file=checkpoint_filename, exist=start_from_checkpoint)
if (start_from_checkpoint) then
    call get_checkpoint
    if (initial_observation_index == no_of_equilibration_sweeps + no_of_observations - 1) then
        ! the last checkpoint was taken at the end of the recorded temperature index
        initial_temperature_index = initial_temperature_index + 1
    end if
    initial_observation_index = initial_observation_index + 1 ! restart at next observation
else
    initial_temperature_index = 0
    initial_observation_index = 0
end if
temperature = temperature + initial_temperature_index * magnitude_of_temperature_increments

! message informing the start of the Markov process
if (algorithm_name == '3dxy-gaussian-noise-metropolis') then
    write(6, '(A, A, A, I4.4, A, I4.4, A, I4.4, A)') 'Starting the ', trim(algorithm_name), ' simulation on ', &
                                                        integer_lattice_length, 'x', integer_lattice_length, 'x', &
                                                        integer_lattice_length, ' lattice sites.'
else
    write(6, '(A, A, A, I4.4, A, I4.4, A)') 'Starting the ', trim(algorithm_name), ' simulation on ', &
                                                integer_lattice_length, 'x', integer_lattice_length, ' lattice sites.'
end if

return
end subroutine pre_simulation_processes
