subroutine output_metropolis_acceptance_rates(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index
double precision :: acceptance_rate_of_field_rotations, acceptance_rate_of_external_global_moves

acceptance_rate_of_field_rotations = no_of_accepted_field_rotations_per_site / dfloat(no_of_samples)

write(filename, '(A, "/temp_", I2.2, "/acceptance_rates.csv")') trim(output_directory), temperature_index
open(unit=20, file = filename)
write(20, '("#", A30, "; ", I0.4, "x", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
        trim(algorithm_name), integer_lattice_length, integer_lattice_length, integer_lattice_length, temperature

write(20, 100) width_of_proposal_interval, acceptance_rate_of_field_rotations
close(20)

100 format(ES24.14, ", ", ES24.14)

return
end subroutine output_metropolis_acceptance_rates
