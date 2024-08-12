subroutine create_sample_files(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index

if (measure_magnetisation) then
    write(filename, '(A, "/temp_", I2.2, "/magnetisation.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 30, temperature_index)
    else
        open(unit=30, file=filename)
        call print_file_header(30)
    end if
end if

if (measure_helicity) then
    write(filename, '(A, "/temp_", I2.2, "/1st_deriv_of_potential.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 40, temperature_index)
    else
        open(unit=40, file=filename)
        call print_file_header(40)
    end if
    write(filename, '(A, "/temp_", I2.2, "/2nd_deriv_of_potential.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 50, temperature_index)
    else
        open(unit=50, file=filename)
        call print_file_header(50)
    end if
end if

if (measure_potential) then
    write(filename, '(A, "/temp_", I2.2, "/potential.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 60, temperature_index)
    else
        open(unit=60, file=filename)
        call print_file_header(60)
    end if
end if

if (measure_potential_minimising_twists) then
    write(filename, '(A, "/temp_", I2.2, "/potential_minimising_twists.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 70, temperature_index)
    else
        open(unit=70, file=filename)
        call print_file_header(70)
    end if
end if

if (measure_external_global_moves) then
    write(filename, '(A, "/temp_", I2.2, "/external_global_moves.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 80, temperature_index)
    else
        open(unit=80, file=filename)
        call print_file_header(80)
    end if
end if

return
end subroutine create_sample_files


subroutine print_file_header(file_index)
use variables
implicit none
integer :: file_index

write(file_index, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
        trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature

return
end subroutine print_file_header
