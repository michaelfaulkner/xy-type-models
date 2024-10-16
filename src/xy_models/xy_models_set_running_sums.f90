subroutine set_running_sums
use variables
implicit none
character(100) :: filename

if (measure_magnetisation) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/magnetic_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) magnetic_norm_sum
        read(11, 200) magnetic_norm_squared_sum
        read(11, 200) magnetic_norm_quartic_sum
        close(11)
    else
        magnetic_norm_sum = 0.0d0
        magnetic_norm_squared_sum = 0.0d0
        magnetic_norm_quartic_sum = 0.0d0
    end if
end if

if (measure_helicity) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/helicity_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) inverse_vacuum_perm_sum
        read(11, 200) inverse_vacuum_perm_squared_sum
        read(11, 200) macro_josephson_current_sum(1)
        read(11, 200) macro_josephson_current_sum(2)
        read(11, 200) macro_josephson_current_squared_sum
        read(11, 200) macro_josephson_current_quartic_sum
        close(11)
    else
        inverse_vacuum_perm_sum = 0.0d0
        inverse_vacuum_perm_squared_sum = 0.0d0
        macro_josephson_current_sum = (/ 0.0d0, 0.0d0 /)
        macro_josephson_current_squared_sum = 0.0d0
        macro_josephson_current_quartic_sum = 0.0d0
    end if
end if

if (measure_potential) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/potential_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) potential_sum
        read(11, 200) potential_squared_sum
        read(11, 200) potential_quartic_sum
        close(11)
    else
        potential_sum = 0.0d0
        potential_squared_sum = 0.0d0
        potential_quartic_sum = 0.0d0
    end if
end if

if (measure_hot_twist_relaxations) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/hot_twist_relaxations_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 100) hot_twist_relaxations_sum(1)
        read(11, 100) hot_twist_relaxations_sum(2)
        read(11, 100) hot_twist_relaxations_squared_sum
        read(11, 100) hot_twist_relaxations_quartic_sum
        close(11)
    else
        hot_twist_relaxations_sum = (/ 0, 0 /)
        hot_twist_relaxations_squared_sum = 0
        hot_twist_relaxations_quartic_sum = 0
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
        read(11, 200) emergent_field_zero_mode_sum(1)
        read(11, 200) emergent_field_zero_mode_sum(2)
        read(11, 200) emergent_field_zero_mode_squared_sum
        read(11, 200) emergent_field_zero_mode_quartic_sum
        read(11, 100) topological_sector_sum(1)
        read(11, 100) topological_sector_sum(2)
        read(11, 100) topological_sector_squared_sum
        read(11, 100) topological_sector_quartic_sum
        close(11)
    else
        emergent_field_zero_mode_sum = (/ 0.0d0, 0.0d0 /)
        emergent_field_zero_mode_squared_sum = 0.0d0
        emergent_field_zero_mode_quartic_sum = 0.0d0
        topological_sector_sum = (/ 0, 0 /)
        topological_sector_squared_sum = 0
        topological_sector_quartic_sum = 0
    end if
end if

100 format(I12)
200 format(ES25.14)

return
end subroutine set_running_sums
