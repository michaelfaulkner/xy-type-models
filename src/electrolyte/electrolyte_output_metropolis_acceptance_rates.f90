subroutine output_metropolis_acceptance_rates
use variables
implicit none
double precision :: acceptance_rate_of_field_rotations, acceptance_rate_of_charge_hops
double precision :: acceptance_rate_of_external_global_moves
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename

acceptance_rate_of_field_rotations = float(no_of_accepted_field_rotations) / (measurements * volume)
acceptance_rate_of_charge_hops = float(no_of_accepted_charge_hops) / (measurements * volume)
acceptance_rate_of_external_global_moves = float(no_of_accepted_external_global_moves) / (measurements * volume)

write(filename, '(A, F4.2, "//acceptance_rates.csv")') trim(output_directory)//trim(temperature_string), temperature
open(unit=30, file = filename)

if (.not.(use_external_global_moves)) then
    write(30, 100) width_of_proposal_interval, acceptance_rate_of_field_rotations, acceptance_rate_of_charge_hops
else
    write(30, 200) width_of_proposal_interval, acceptance_rate_of_field_rotations, acceptance_rate_of_charge_hops, &
                        acceptance_rate_of_external_global_moves
end if
close(30)

100 format(ES24.14, ", ", ES24.14, ", ", ES24.14)
200 format(ES24.14, ", ", ES24.14, ", ", ES24.14, ", ", ES24.14)

return
end subroutine output_metropolis_acceptance_rates
