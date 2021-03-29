! **************************************
! Global twist (Metropolis) sampling
! **************************************

subroutine global_twist_XY
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

     do i = 1, sites
        deltaE = deltaE - cos(theta(i) - theta(neg_x(i)) + rand1 * twopi / side) + cos(theta(i) - theta(neg_x(i)))
     end do

     if ((deltaE < 0) .or. (rand() < exp(-beta * deltaE))) then
        rand2 = int(side * rand())                                               ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
        do i = 1, sites
           theta(i) = mod(theta(i) + rand1 * twopi * mod(i + rand2, side) / side, twopi)
        end do
        accept_twist = accept_twist + 1
     end if

  else                                                                               ! Y TWIST

     rand1 = floor(2 * rand())                                                       ! CHOOSES +VE OR -VE TWIST DIRECTION
     if (rand1 == 0) then
        rand1 = -1
     end if     

     do i = 1, sites
        deltaE = deltaE - cos(theta(i) - theta(neg_y(i)) + rand1 * twopi / side) + cos(theta(i) - theta(neg_y(i)))
     end do

     if ((deltaE < 0) .or. (rand() < exp(-beta * deltaE))) then
        rand2 = int(side * rand())                                                   ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
        do i = 1, sites
           theta(i) = mod(theta(i) + rand1 * twopi * (int(i / side) + rand2) / side, twopi)
        end do
        accept_twist = accept_twist + 1
     end if

  end if

end subroutine global_twist_XY
