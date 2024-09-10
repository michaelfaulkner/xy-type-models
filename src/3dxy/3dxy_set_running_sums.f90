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

100 format(I12)
200 format(ES25.14)

return
end subroutine set_running_sums
