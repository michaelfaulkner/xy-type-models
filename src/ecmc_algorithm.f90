program ecmc_algorithm
use variables
implicit none
integer :: i, j
double precision :: start_time, end_time, no_of_events_per_sweep

call pre_simulation_processes
call cpu_time(start_time)

do i = 0, no_of_temperature_increments
    write(6, '(A, F4.2)') 'Temperature = ', temperature
    beta = 1.0d0 / temperature
    spin_space_distance_between_observations = dfloat(no_of_sites)
    call create_sample_file

    do j = 1, no_of_equilibration_sweeps
        call single_event_chain
        ! spin_space_distance_between_observations adaptor
        if (mod(j, no_of_equilibration_sweeps / 10) == 0) then
            no_of_events_per_sweep = 1.0d-2 * dfloat(no_of_events)
            if (no_of_events_per_sweep > 1.1d0 * dfloat(no_of_sites)) then
                spin_space_distance_between_observations = 9.99000999000999d-1 * spin_space_distance_between_observations
            else if (no_of_events_per_sweep < 9.09090909090909d-1 * dfloat(no_of_sites)) then
                spin_space_distance_between_observations = 1.001d0 * spin_space_distance_between_observations
            end if
        end if
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call get_and_print_observation
    end do

    no_of_events = 0
    no_of_accepted_external_global_moves = 0

    do j = 1, no_of_observations
        call single_event_chain
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call get_and_print_observation
    end do

    call output_no_of_events
    temperature = temperature + magnitude_of_temperature_increments
    call initialise_field_configuration(.false.)
end do

call cpu_time(end_time)
write(6, '(A, ES9.3, A)') 'Markov process complete.  Total runtime of the simulation = ', end_time - start_time, &
                            ' seconds.'

end program ecmc_algorithm


subroutine output_no_of_events
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename

write (filename, '(A, F4.2, "//no_of_events.csv")' ) trim(output_directory)//trim(temperature_string), temperature
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
