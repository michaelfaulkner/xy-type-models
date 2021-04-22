subroutine attempt_external_global_move
  use variables
  implicit none
  integer i, plusMinus
  double precision EsumNew, deltaU

  call calculate_electric_field_sum
  deltaU = 0
  
  if (floor(2 * rand()) .eq. 0) then
     
     plusMinus = 2 * floor(2 * rand()) - 1
     EsumNew = Esum_x + plusMinus * elementaryCharge * side
     deltaU = 0.5 * (EsumNew * EsumNew - Esum_x * Esum_x) / volume

     ! METROPOLIS FILTER

     if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
        do i = 1, sites
           electric_field_x(i) = electric_field_x(i) + plusMinus * elementaryCharge / side
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
     end if
     
  else
     
     plusMinus = 2 * floor(2 * rand()) - 1
     EsumNew = Esum_y + plusMinus * elementaryCharge * side
     deltaU = 0.5 * (EsumNew * EsumNew - Esum_y * Esum_y) / volume
          
     ! METROPOLIS FILTER

     if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
        do i = 1, sites
           electric_field_y(i) = electric_field_y(i) + plusMinus * elementaryCharge / side
        end do
        no_of_accepted_external_global_moves = no_of_accepted_external_global_moves + 1
     end if
     
  end if
  
end subroutine attempt_external_global_move
