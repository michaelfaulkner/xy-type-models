program metropolis_algorithm
use variables
implicit none
integer :: temperature_index, observation_index
real :: start_time, end_time
double precision :: acceptance_rate_of_field_rotns

call pre_simulation_processes
if (.not.start_from_checkpoint) then
    call reset_metropolis_acceptance_counters
end if

if (.not.simulation_complete) then
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
                call metropolis_sweep
                ! step-size adaptor
                if (mod(observation_index, 100) == 0) then
                    acceptance_rate_of_field_rotns = 1.0d-2 * no_of_accepted_field_rotations_per_site / &
                                                            (1.0d0 - charge_hop_proportion)
                    if (acceptance_rate_of_field_rotns > target_acceptance_rate_of_field_rotations + 5.0d-2) then
                        if (width_of_proposal_interval < 10.0d0) then ! to avoid huge step sizes at high temperature
                            width_of_proposal_interval = 1.1d0 * width_of_proposal_interval
                        end if
                    else if (acceptance_rate_of_field_rotns < target_acceptance_rate_of_field_rotations - 5.0d-2) then
                        width_of_proposal_interval = 9.0d-1 * width_of_proposal_interval
                    end if
                    no_of_accepted_field_rotations_per_site = 0.0d0
                end if
                call get_and_print_observation
                call do_checkpointing(temperature_index, observation_index)
            end do
            call reset_metropolis_acceptance_counters
            initial_observation_index = no_of_equilibration_sweeps
        end if

        do observation_index = initial_observation_index, no_of_equilibration_sweeps + no_of_observations - 1
            if (use_external_global_moves) then
                call attempt_external_global_moves
            end if
            call metropolis_sweep
            call get_and_print_observation
            if (observation_index < no_of_equilibration_sweeps + no_of_observations - 1) then
                call do_checkpointing(temperature_index, observation_index) ! we don't checkpoint at end of temp index
            end if
        end do

        call output_metropolis_acceptance_rates(temperature_index)
        call reset_metropolis_acceptance_counters
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

end program metropolis_algorithm
