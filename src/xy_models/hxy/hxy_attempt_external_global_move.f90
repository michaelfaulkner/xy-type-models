! **************************************
! Global twist (Metropolis) sampling
! **************************************

subroutine attempt_external_global_move
  use variables
  implicit none
  integer i, rand1, rand2
  double precision deltaE, diff

  call calculate_emergent_field
  deltaE = 0
  
  if (floor(2 * rand()) .eq. 0) then                                                 ! Ebar_x MOVE
     
     rand1 = floor(2 * rand())                                                       ! CHOOSES +VE OR -VE TWIST DIRECTION
     if (rand1 .eq. 0) then
        rand1 = -1
     end if
     
     do i = 1, no_of_sites
        diff = spin_field(i) - spin_field(neg_y(i)) + rand1 * twopi / integer_lattice_length
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        deltaE = deltaE + 0.5 * diff ** 2 - 0.5 * emergent_field_x(i) ** 2
     end do
     
     if ((deltaE < 0.0d0) .or. (exp(-beta * deltaE) .gt. rand())) then
        rand2 = int(integer_lattice_length * rand())                                               ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
        do i = 1, no_of_sites
           spin_field(i) = mod(spin_field(i) + rand1 * twopi * (int(i / integer_lattice_length) + rand2) / &
                                 integer_lattice_length, twopi)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
     end if
     
  else                                                                            ! Ebar_Y MOVE
     
     rand1 = floor(2 * rand())                                                    ! CHOOSES +VE OR -VE TWIST DIRECTION
     if (rand1 .eq. 0) then
        rand1 = -1
     end if

     do i = 1, no_of_sites
        diff = -(spin_field(i) - spin_field(neg_x(i))) + rand1 * twopi / integer_lattice_length
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        deltaE = deltaE + 0.5 * diff ** 2 - 0.5 * emergent_field_y(i) ** 2
     end do
     
     if ((deltaE < 0.0d0) .or. (exp(-beta * deltaE) .gt. rand())) then
        rand2 = int(integer_lattice_length * rand())                                               ! TO PICK A RANDOM STARTING LINE FROM WHICH TO TWIST
        do i = 1, no_of_sites
           spin_field(i) = mod(spin_field(i) - rand1 * twopi * mod(i + rand2, integer_lattice_length) / &
                                 integer_lattice_length, twopi)
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
     end if
     
  end if
  
end subroutine attempt_external_global_move
