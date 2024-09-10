subroutine record_running_sums
use variables
implicit none
character(100) :: temporary_filename, filename

if (measure_magnetisation) then
    write(filename, '(A, "/magnetic_running_sums.csv")') trim(output_directory)
    write(temporary_filename, '(A, "/temporary_running_sums.csv")') trim(output_directory)
    open(unit=11, file=temporary_filename)
    write(11, 200) magnetic_norm_sum
    write(11, 200) magnetic_norm_squared_sum
    write(11, 200) magnetic_norm_quartic_sum
    close(11)
    call system('mv ' // trim(temporary_filename) // ' ' // trim(filename))
end if

if (measure_helicity) then
    write(filename, '(A, "/helicity_running_sums.csv")') trim(output_directory)
    write(temporary_filename, '(A, "/temporary_running_sums.csv")') trim(output_directory)
    open(unit=11, file=temporary_filename)
    write(11, 200) inverse_vacuum_perm_sum
    write(11, 200) inverse_vacuum_perm_squared_sum
    write(11, 200) macro_josephson_current_sum(1)
    write(11, 200) macro_josephson_current_sum(2)
    write(11, 200) macro_josephson_current_squared_sum
    write(11, 200) macro_josephson_current_quartic_sum
    close(11)
    call system('mv ' // trim(temporary_filename) // ' ' // trim(filename))
end if

if (measure_potential) then
    write(filename, '(A, "/potential_running_sums.csv")') trim(output_directory)
    open(unit=11, file=filename)
    write(11, 200) potential_sum
    write(11, 200) potential_squared_sum
    write(11, 200) potential_quartic_sum
    close(11)
end if

if (measure_hot_twist_relaxations) then
    write(filename, '(A, "/hot_twist_relaxations_running_sums.csv")') trim(output_directory)
    write(temporary_filename, '(A, "/temporary_running_sums.csv")') trim(output_directory)
    open(unit=11, file=temporary_filename)
    write(11, 100) hot_twist_relaxations_sum(1)
    write(11, 100) hot_twist_relaxations_sum(2)
    write(11, 100) hot_twist_relaxations_squared_sum
    write(11, 100) hot_twist_relaxations_quartic_sum
    close(11)
    call system('mv ' // trim(temporary_filename) // ' ' // trim(filename))
end if

if (measure_twist_relaxations) then
    write(filename, '(A, "/twist_relaxations_running_sums.csv")') trim(output_directory)
    write(temporary_filename, '(A, "/temporary_running_sums.csv")') trim(output_directory)
    open(unit=11, file=temporary_filename)
    write(11, 100) twist_relaxations_sum(1)
    write(11, 100) twist_relaxations_sum(2)
    write(11, 100) twist_relaxations_squared_sum
    write(11, 100) twist_relaxations_quartic_sum
    close(11)
    call system('mv ' // trim(temporary_filename) // ' ' // trim(filename))
end if

if (measure_emergent_field) then
    write(filename, '(A, "/emergent_field_running_sums.csv")') trim(output_directory)
    write(temporary_filename, '(A, "/temporary_running_sums.csv")') trim(output_directory)
    open(unit=11, file=temporary_filename)
    write(11, 200) emergent_field_zero_mode_sum(1)
    write(11, 200) emergent_field_zero_mode_sum(2)
    write(11, 200) emergent_field_zero_mode_squared_sum
    write(11, 200) emergent_field_zero_mode_quartic_sum
    write(11, 100) topological_sector_sum(1)
    write(11, 100) topological_sector_sum(2)
    write(11, 100) topological_sector_squared_sum
    write(11, 100) topological_sector_quartic_sum
    close(11)
    call system('mv ' // trim(temporary_filename) // ' ' // trim(filename))
end if

100 format(I12)
200 format(ES25.14)

return
end subroutine record_running_sums
