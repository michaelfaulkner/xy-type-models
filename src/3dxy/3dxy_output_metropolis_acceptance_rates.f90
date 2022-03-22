subroutine output_metropolis_acceptance_rates(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index
double precision :: acceptance_rate_of_field_rotations, acceptance_rate_of_external_global_moves

acceptance_rate_of_field_rotations = dfloat(no_of_accepted_field_rotations) / dfloat(no_of_observations) &
                                        / dfloat(no_of_sites)

write(filename, '(A, "/temp_", I2.2, "/acceptance_rates.csv")') trim(output_directory), temperature_index
open(unit=30, file = filename)

write(30, 100) width_of_proposal_interval, acceptance_rate_of_field_rotations
close(30)

100 format(ES24.14, ", ", ES24.14)

return
end subroutine output_metropolis_acceptance_rates
