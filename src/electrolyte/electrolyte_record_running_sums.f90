subroutine record_running_sums
use variables
implicit none
character(100) :: temporary_filename, filename

if (measure_electric_field_sum) then
    write(filename, '(A, "/electric_field_running_sums.csv")') trim(output_directory)
    write(temporary_filename, '(A, "/temporary_running_sums.csv")') trim(output_directory)
    open(unit=11, file=temporary_filename)
    write(11, 200) raw_electric_field_zero_mode_sum(1)
    write(11, 200) raw_electric_field_zero_mode_sum(2)
    write(11, 200) raw_electric_field_zero_mode_squared_sum
    write(11, 200) raw_electric_field_zero_mode_quartic_sum
    write(11, 100) topological_sector_sum(1)
    write(11, 100) topological_sector_sum(2)
    write(11, 100) topological_sector_squared_sum
    write(11, 100) topological_sector_quartic_sum
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

100 format(I12)
200 format(ES25.14)

return
end subroutine record_running_sums
