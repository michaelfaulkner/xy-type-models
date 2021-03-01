! **************************************
! Global twist (Metropolis) sampling
! **************************************

subroutine global_twist_HXY
  use variables
  implicit none
  integer i, rand1, rand2
  real*8 deltaE, diff

  call top_field
  deltaE = 0
  
  if (floor(2 * rand()) .eq. 0) then                                                 ! Ebar_x MOVE
     
     rand1 = floor(2 * rand())                                                       ! CHOOSES +VE OR -VE TWIST DIRECTION
     if (rand1 .eq. 0) then
        rand1 = -1
     end if
     
     do i = 0, sites - 1
        diff = theta(i) - theta(neg_y(i)) + rand1 * twopi / side
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        deltaE = deltaE + 0.5 * diff ** 2 - 0.5 * top_x(i) ** 2
     end do
     
     if ((deltaE .lt. -epsilon) .or. (exp(-beta * deltaE) .gt. rand())) then
        rand2 = int(side * rand())                                               ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
        do i = 0, sites - 1
           theta(i) = mod(theta(i) + rand1 * twopi * (int(i / side) + rand2) / side, twopi)
        end do
        accept_twist = accept_twist + 1
     end if
     
  else                                                                            ! Ebar_Y MOVE
     
     rand1 = floor(2 * rand())                                                    ! CHOOSES +VE OR -VE TWIST DIRECTION
     if (rand1 .eq. 0) then
        rand1 = -1
     end if

     do i = 0, sites - 1
        diff = -(theta(i) - theta(neg_x(i))) + rand1 * twopi / side
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        deltaE = deltaE + 0.5 * diff ** 2 - 0.5 * top_y(i) ** 2
     end do
     
     if ((deltaE .lt. -epsilon) .or. (exp(-beta * deltaE) .gt. rand())) then
        rand2 = int(side * rand())                                               ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
        do i = 0, sites - 1
           theta(i) = mod(theta(i) - rand1 * twopi * mod(i + rand2, side) / side, twopi)
        end do
        accept_twist = accept_twist + 1
     end if
     
  end if
  
end subroutine global_twist_HXY
