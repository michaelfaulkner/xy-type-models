program ecmc_algorithm
use variables
implicit none
character(100) :: config_file
integer :: i, j, seed, start
double precision :: magnitude_of_temperature_increments

! verify that the something has been parsed to the exectuable
if (command_argument_count() /= 1) then
    write(6, *) 'Error: parse configuration file to executable'
    stop
end if
! read in config file
call get_command_argument(1, config_file)
open (unit=1, file=config_file)

call input(seed, start)
call setup_periodic_boundaries
call create_sample_files
call initialise_field_configuration(start)
call randinit(seed)
write(6, '(A, F16.14)') 'Initial random number = ', rand(seed)

if (no_of_temperature_increments == 0) then
    magnitude_of_temperature_increments = 0.0
else
    magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
end if

do i = 0, no_of_temperature_increments
    write(6, '(A, ES8.2)') 'Temperature = ', temperature
    beta = 1.0 / temperature

    do j = 1, therm_sweeps
        call single_event_chain
        if (twist == 1) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    no_of_events = 0
    no_of_accepted_external_global_moves = 0

    do j = 1, measurements
        call single_event_chain
        if (twist == 1) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    call output_no_of_events
    temperature = temperature + magnitude_of_temperature_increments
end do
end program ecmc_algorithm


subroutine output_no_of_events
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename

write (filename, '(A, F4.2, "//number_of_events.dat")' ) trim(output_directory)//trim(temperature_string), temperature
open(unit = 300, file = filename)
write(300, 100) no_of_events
if (twist == 1) then
    write(300, 100) no_of_accepted_external_global_moves
end if
close(300)

100 format(I20)

return
end subroutine output_no_of_events
