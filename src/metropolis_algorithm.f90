program metropolis_algorithm
use variables
implicit none
integer :: temperature_index, observation_index
double precision :: start_time, end_time, acceptance_rate_of_field_rotations

call pre_simulation_processes
call cpu_time(start_time)

do temperature_index = 0, no_of_temperature_increments
    write(6, '(A, F6.4)') 'Temperature = ', temperature
    beta = 1.0d0 / temperature
    call create_sample_files(temperature_index)

    call reset_metropolis_acceptance_counters
    do observation_index = 1, no_of_equilibration_sweeps
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call metropolis_sweep
        ! step-size adaptor
        if (mod(observation_index, 100) == 0) then
            acceptance_rate_of_field_rotations = 1.0d-2 * dfloat(no_of_accepted_field_rotations) / dfloat(no_of_sites) &
                                                    / (1.0d0 - charge_hop_proportion)
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
    end do

    call reset_metropolis_acceptance_counters
    do observation_index = 1, no_of_observations
        if (use_external_global_moves) then
            call attempt_external_global_moves
        end if
        call metropolis_sweep
        call get_and_print_observation
    end do

    call output_metropolis_acceptance_rates(temperature_index)
    temperature = temperature + magnitude_of_temperature_increments
    call initialise_field_configuration(.false.)
end do

call cpu_time(end_time)
write(6, '(A, ES9.3, A)') 'Markov process complete.  Total runtime of the simulation = ', end_time - start_time, &
                            ' seconds.'

end program metropolis_algorithm
