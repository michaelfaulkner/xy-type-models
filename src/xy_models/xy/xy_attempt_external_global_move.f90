subroutine attempt_external_global_move
use variables
implicit none
integer :: i, rand1, rand2
double precision :: deltaE

deltaE = 0

if (floor(2 * rand()) == 0) then                                                 ! X TWIST

 rand1 = floor(2 * rand())                                                       ! CHOOSES +VE OR -VE TWIST DIRECTION
 if (rand1 == 0) then
    rand1 = -1
 end if

 do i = 1, no_of_sites
    deltaE = deltaE - cos(theta(i) - theta(neg_x(i)) + rand1 * twopi / integer_lattice_length) + cos(theta(i) - theta(neg_x(i)))
 end do

 if ((deltaE < 0) .or. (rand() < exp(-beta * deltaE))) then
    rand2 = int(integer_lattice_length * rand())                                               ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
    do i = 1, no_of_sites
       theta(i) = mod(theta(i) + rand1 * twopi * mod(i + rand2, integer_lattice_length) / integer_lattice_length, twopi)
    end do
    no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
 end if

else                                                                               ! Y TWIST

 rand1 = floor(2 * rand())                                                       ! CHOOSES +VE OR -VE TWIST DIRECTION
 if (rand1 == 0) then
    rand1 = -1
 end if

 do i = 1, no_of_sites
    deltaE = deltaE - cos(theta(i) - theta(neg_y(i)) + rand1 * twopi / integer_lattice_length) + cos(theta(i) - theta(neg_y(i)))
 end do

 if ((deltaE < 0) .or. (rand() < exp(-beta * deltaE))) then
    rand2 = int(integer_lattice_length * rand())                                                   ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
    do i = 1, no_of_sites
       theta(i) = mod(theta(i) + rand1 * twopi * (int(i / integer_lattice_length) + rand2) / integer_lattice_length, twopi)
    end do
    no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
 end if

end if

end subroutine attempt_external_global_move
