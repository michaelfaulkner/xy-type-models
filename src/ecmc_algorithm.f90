program ecmc_algorithm
use variables
implicit none
integer :: temperature_index, observation_index
double precision :: start_time, end_time, no_of_events_per_sweep

call pre_simulation_processes
call cpu_time(start_time)

do temperature_index = 0, no_of_temperature_increments
    write(6, '(A, F6.4)') 'Temperature = ', temperature
    beta = 1.0d0 / temperature
    spin_space_distance_between_observations = 1.0d0 * dfloat(no_of_sites)
    call create_sample_files(temperature_index)

    do observation_index = 1, no_of_equilibration_sweeps
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call single_event_chain
        call get_and_print_observation
    end do

    no_of_events = 0
    no_of_accepted_external_global_moves = 0

    do observation_index = 1, no_of_observations
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call single_event_chain
        call get_and_print_observation
    end do

    call output_no_of_events(temperature_index)
    temperature = temperature + magnitude_of_temperature_increments
    call initialise_field_configuration(.false.)
end do

call cpu_time(end_time)
write(6, '(A, ES9.3, A)') 'Markov process complete.  Total runtime of the simulation = ', end_time - start_time, &
                            ' seconds.'

end program ecmc_algorithm


subroutine output_no_of_events(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index

write (filename, '(A, "/temp_", I2.2, "/no_of_events.csv")' ) trim(output_directory), temperature_index
open(unit=300, file = filename)
if (.not.(use_external_global_moves)) then
    write(300, 100) no_of_events
else
    write(300, 200) no_of_events, 0.5d0 * dfloat(no_of_accepted_external_global_moves) / dfloat(no_of_observations)
end if
close(300)

100 format(I20)
200 format(I20, ", ", ES24.14)

return
end subroutine output_no_of_events
