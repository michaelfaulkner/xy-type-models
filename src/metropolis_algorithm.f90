program metropolis_algorithm
use variables
implicit none
integer :: temperature_index, observation_index
real :: start_time, end_time
double precision :: acceptance_rate_of_field_rotations

call pre_simulation_processes
call cpu_time(start_time)

do temperature_index = initial_temperature_index, no_of_temperature_increments
    write(6, '(A, F6.4)') 'Temperature = ', temperature
    beta = 1.0d0 / temperature
    call create_sample_files(temperature_index)

    if (initial_observation_index < no_of_equilibration_sweeps) then
        call reset_metropolis_acceptance_counters
        do observation_index = initial_observation_index, no_of_equilibration_sweeps - 1
            if (use_external_global_moves) then
                call attempt_external_global_moves
            end if
            call metropolis_sweep
            ! step-size adaptor
            if (mod(observation_index, 100) == 0) then
                acceptance_rate_of_field_rotations = 1.0d-2 * dfloat(no_of_accepted_field_rotations) / &
                                                        dfloat(no_of_sites) / (1.0d0 - charge_hop_proportion)
                if (acceptance_rate_of_field_rotations > target_acceptance_rate_of_field_rotations + 5.0d-2) then
                    if (width_of_proposal_interval < 10.0d0) then ! to avoid huge step sizes at high temperature
                        width_of_proposal_interval = 1.1d0 * width_of_proposal_interval
                    end if
                else if (acceptance_rate_of_field_rotations < target_acceptance_rate_of_field_rotations - 5.0d-2) then
                    width_of_proposal_interval = 9.0d-1 * width_of_proposal_interval
                end if
                no_of_accepted_field_rotations = 0
            end if
            call get_and_print_observation
            call do_checkpointing(temperature_index, observation_index)
        end do
        initial_observation_index = no_of_equilibration_sweeps
    end if

    call reset_metropolis_acceptance_counters
    do observation_index = initial_observation_index, no_of_equilibration_sweeps + no_of_observations - 1
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call metropolis_sweep
        call get_and_print_observation
        call do_checkpointing(temperature_index, observation_index)
    end do

    call output_metropolis_acceptance_rates(temperature_index)
    start_from_checkpoint = .false.
    initial_observation_index = 0
    temperature = temperature + magnitude_of_temperature_increments
    call initialise_field_configuration(.false.)
end do

call delete_last_checkpoint
call cpu_time(end_time)
write(6, '(A, ES9.3, A)') 'Markov process complete.  Total runtime of the simulation = ', &
                            (end_time - start_time) / 60.0, ' minutes.'

end program metropolis_algorithm
