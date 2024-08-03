program ecmc_algorithm
use variables
implicit none
integer :: temperature_index, sample_index
real :: start_time, end_time

call pre_simulation_processes
if (.not.start_from_checkpoint) then
    call reset_event_counters
end if

if (.not.simulation_complete) then
    spin_space_distance_between_samples = 1.0d0 * dfloat(no_of_sites)
    call cpu_time(start_time)
    do temperature_index = initial_temperature_index, no_of_temperature_increments
        write(6, '(A, F6.4)') 'Temperature = ', temperature
        beta = 1.0d0 / temperature
        if (print_samples) then
            call create_sample_files(temperature_index)
        end if
        call set_running_sums

        if (initial_sample_index < no_of_equilibration_sweeps) then
            do sample_index = initial_sample_index, no_of_equilibration_sweeps - 1
                if (use_external_global_moves) then
                    call attempt_external_global_moves
                end if
                call single_event_chain
                if (print_samples) then
                    call process_sample(sample_index)
                end if
                call do_checkpointing(temperature_index, sample_index)
            end do
            call reset_event_counters
            initial_sample_index = no_of_equilibration_sweeps
        end if

        do sample_index = initial_sample_index, no_of_equilibration_sweeps + no_of_samples - 1
            if (use_external_global_moves) then
                call attempt_external_global_moves
            end if
            call single_event_chain
            call process_sample(sample_index)
            if (sample_index < no_of_equilibration_sweeps + no_of_samples - 1) then
                call do_checkpointing(temperature_index, sample_index) ! we don't checkpoint at end of temp index
            end if
        end do

        call output_event_rate(temperature_index)
        call reset_event_counters
        call output_summary_statistics(temperature_index)
        start_from_checkpoint = .false.
        initial_sample_index = 0
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
    write(300, 100) no_of_events_per_unit_spin_space_distance / dfloat(no_of_samples)
else
    write(300, 200) no_of_events_per_unit_spin_space_distance / dfloat(no_of_samples), &
                        0.5d0 * no_of_accepted_external_global_moves / dfloat(no_of_samples)
end if
close(300)

100 format(ES24.14)
200 format(ES24.14, ", ", ES24.14)

return
end subroutine output_event_rate
