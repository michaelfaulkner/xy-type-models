subroutine set_running_sums
use variables
implicit none
character(100) :: filename

if (measure_magnetisation) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/magnetic_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) raw_magnetic_norm_sum
        read(11, 200) raw_magnetic_norm_squared_sum
        read(11, 200) raw_magnetic_norm_quartic_sum
        close(11)
    else
        raw_magnetic_norm_sum = 0.0d0
        raw_magnetic_norm_squared_sum = 0.0d0
        raw_magnetic_norm_quartic_sum = 0.0d0
    end if
end if

if (measure_helicity) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/helicity_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) raw_inverse_vacuum_perm_sum
        read(11, 200) raw_inverse_vacuum_perm_squared_sum
        read(11, 200) raw_macro_josephson_current_sum(1)
        read(11, 200) raw_macro_josephson_current_sum(2)
        read(11, 200) raw_macro_josephson_current_squared_sum
        read(11, 200) raw_macro_josephson_current_quartic_sum
        close(11)
    else
        raw_inverse_vacuum_perm_sum = 0.0d0
        raw_inverse_vacuum_perm_squared_sum = 0.0d0
        raw_macro_josephson_current_sum = (/ 0.0d0, 0.0d0 /)
        raw_macro_josephson_current_squared_sum = 0.0d0
        raw_macro_josephson_current_quartic_sum = 0.0d0
    end if
end if

if (measure_potential) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/potential_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) potential_sum
        read(11, 200) potential_squared_sum
        close(11)
    else
        potential_sum = 0.0d0
        potential_squared_sum = 0.0d0
    end if
end if

if (measure_potential_minimising_twists) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/potential_minimising_twists_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 100) potential_minimising_twists_sum(1)
        read(11, 100) potential_minimising_twists_sum(2)
        read(11, 100) potential_minimising_twists_squared_sum
        read(11, 100) potential_minimising_twists_quartic_sum
        close(11)
    else
        potential_minimising_twists_sum = (/ 0, 0 /)
        potential_minimising_twists_squared_sum = 0
        potential_minimising_twists_quartic_sum = 0
    end if
end if

if (measure_twist_relaxations) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/twist_relaxations_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 100) twist_relaxations_sum(1)
        read(11, 100) twist_relaxations_sum(2)
        read(11, 100) twist_relaxations_squared_sum
        read(11, 100) twist_relaxations_quartic_sum
        close(11)
    else
        twist_relaxations_sum = (/ 0, 0 /)
        twist_relaxations_squared_sum = 0
        twist_relaxations_quartic_sum = 0
    end if
end if

if (measure_emergent_field) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/emergent_field_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) sum_of_emergent_field_sum(1)
        read(11, 200) sum_of_emergent_field_sum(2)
        read(11, 200) sum_of_emergent_field_squared_sum
        read(11, 200) sum_of_emergent_field_quartic_sum
        read(11, 100) topological_sector_sum(1)
        read(11, 100) topological_sector_sum(2)
        read(11, 100) topological_sector_squared_sum
        read(11, 100) topological_sector_quartic_sum
        close(11)
    else
        sum_of_emergent_field_sum = (/ 0.0d0, 0.0d0 /)
        sum_of_emergent_field_squared_sum = 0.0d0
        sum_of_emergent_field_quartic_sum = 0.0d0
        topological_sector_sum = (/ 0, 0 /)
        topological_sector_squared_sum = 0
        topological_sector_quartic_sum = 0
    end if
end if

100 format(I12)
200 format(ES25.14)

return
end subroutine set_running_sums
