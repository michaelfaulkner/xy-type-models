program ecmc_algorithm
use variables
implicit none
integer :: temperature_index, observation_index
real :: start_time, end_time

call pre_simulation_processes
if (.not.start_from_checkpoint) then
    call reset_event_counters
end if

if (.not.simulation_complete) then
    spin_space_distance_between_observations = 1.0d0 * dfloat(no_of_sites)
    call cpu_time(start_time)
    do temperature_index = initial_temperature_index, no_of_temperature_increments
        write(6, '(A, F6.4)') 'Temperature = ', temperature
        beta = 1.0d0 / temperature
        call create_sample_files(temperature_index)

        if (initial_observation_index < no_of_equilibration_sweeps) then
            do observation_index = initial_observation_index, no_of_equilibration_sweeps - 1
                if (use_external_global_moves) then
                    call attempt_external_global_moves
                end if
                call single_event_chain
                call get_and_print_observation
                call do_checkpointing(temperature_index, observation_index)
            end do
            call reset_event_counters
            initial_observation_index = no_of_equilibration_sweeps
        end if

        do observation_index = initial_observation_index, no_of_equilibration_sweeps + no_of_observations - 1
            if (use_external_global_moves) then
                call attempt_external_global_moves
            end if
            call single_event_chain
            call get_and_print_observation
            call do_checkpointing(temperature_index, observation_index)
        end do

        call output_event_rate(temperature_index)
        call reset_event_counters
        start_from_checkpoint = .false.
        initial_observation_index = 0
        temperature = temperature + magnitude_of_temperature_increments
        call initialise_field_configuration(.false.)
    end do

    call delete_last_checkpoint
    call cpu_time(end_time)
    write(6, '(A, ES9.3, A)') 'Markov process complete.  Total runtime of the simulation = ', &
                                (end_time - start_time) / 60.0, ' minutes.'
end if

end program ecmc_algorithm


subroutine output_event_rate(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index

write (filename, '(A, "/temp_", I2.2, "/event_rate.csv")' ) trim(output_directory), temperature_index
open(unit=300, file = filename)
if (.not.(use_external_global_moves)) then
    write(300, 100) no_of_events_per_unit_spin_space_distance / dfloat(no_of_observations)
else
    write(300, 200) no_of_events_per_unit_spin_space_distance / dfloat(no_of_observations), &
                        0.5d0 * no_of_accepted_external_global_moves / dfloat(no_of_observations)
end if
close(300)

100 format(I20)
200 format(I20, ", ", ES24.14)

return
end subroutine output_event_rate
