! **************************************
! Global TSF sampling
! **************************************

subroutine markov_chain_TSF_GLE
  use variables
  implicit none
  integer i, plusMinus
  real*8 EsumNew, deltaU

  call measure_Esum
  deltaU = 0
  
  if (floor(2 * rand()) .eq. 0) then
     
     plusMinus = 2 * floor(2 * rand()) - 1
     EsumNew = Esum_x + plusMinus * elementaryCharge * side
     deltaU = 0.5 * (EsumNew * EsumNew - Esum_x * Esum_x) / volume

     ! METROPOLIS FILTER

     if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
        do i = 1, sites
           Efield_x(i) = Efield_x(i) + plusMinus * elementaryCharge / side
        end do
        accept_TSF = accept_TSF + 1
     end if
     
  else
     
     plusMinus = 2 * floor(2 * rand()) - 1
     EsumNew = Esum_y + plusMinus * elementaryCharge * side
     deltaU = 0.5 * (EsumNew * EsumNew - Esum_y * Esum_y) / volume
          
     ! METROPOLIS FILTER

     if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
        do i = 1, sites
           Efield_y(i) = Efield_y(i) + plusMinus * elementaryCharge / side
        end do
        accept_TSF = accept_TSF + 1
     end if
     
  end if
  
end subroutine markov_chain_TSF_GLE
