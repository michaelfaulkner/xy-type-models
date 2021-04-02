subroutine setup_periodic_boundaries
use variables
implicit none
integer :: i

! mod(i + side + 1, side) rather than mod(i + 1, side) below as mod(i, side) doesn't return value in [0, side) if i < 0
do i = 1, sites
    pos_x(i) = i + mod(i + side + 1, side) - mod(i + side, side)
    neg_x(i) = i + mod(i + side - 1, side) - mod(i + side, side)
    pos_y(i) = i + (mod(int(i / side) + side + 1, side) - mod(int(i / side) + side, side)) * side
    neg_y(i) = i + (mod(int(i / side) + side - 1, side) - mod(int(i / side) + side, side)) * side
end do

return
end subroutine setup_periodic_boundaries

subroutine initialise_spin_configuration(start)
use variables
implicit none
integer :: i, start

if (start == 0) then
    do i = 1, sites
        theta(i) = 0.0d0
    end do
else
    do i = 1, sites
        theta(i) = twopi * (rand() - 0.5d0)
    end do
end if

return
end subroutine initialise_spin_configuration

! **************************************
! CREATE NEW DIRECTORY AND FILES FOR NEW TEMP
! **************************************

subroutine initial_measure
use variables
implicit none
character(100) :: filename
character(100) :: temperature_directory
character(100), parameter :: temperature_string="/temp_eq_"

! OPENS NEW DIRECTORY IN WHICH TO SAVE THE MARKOV CHAIN FOR THE CURRENT TEMPERATURE
write (temperature_directory, '(A, F4.2)' ) trim(output_directory)//trim(temperature_string), T
call system ( 'mkdir -p ' // temperature_directory )

write (filename, '(A, F4.2, "//magnetisation_x_sample.dat")' ) trim(output_directory)//trim(temperature_string), &
        T
open(unit=10,file=filename)
write (filename, '(A, F4.2, "//magnetisation_y_sample.dat")' ) trim(output_directory)//trim(temperature_string), &
        T
open(unit=11,file=filename)
write (filename, '(A, F4.2, "//mean_1st_derivative_of_potential_x_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), T
open(unit=12,file=filename)
write (filename, '(A, F4.2, "//mean_1st_derivative_of_potential_y_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), T
open(unit=13,file=filename)
write (filename, '(A, F4.2, "//mean_2nd_derivative_of_potential_x_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), T
open(unit=14,file=filename)
write (filename, '(A, F4.2, "//mean_2nd_derivative_of_potential_y_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), T
open(unit=15,file=filename)
write (filename, '(A, F4.2, "//potential_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=16, file = filename)
write (filename, '(A, F4.2, "//external_minimising_twist_field_x_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), T
open(unit=17, file = filename)
write (filename, '(A, F4.2, "//external_minimising_twist_field_y_sample.dat")' ) &
        trim(output_directory)//trim(temperature_string), T
open(unit=18, file = filename)

return
end subroutine initial_measure
