program test
  integer i

  do i = 0, 20
     write(6,*) 2 * floor(2 * rand()) - 1
  end do

end program test
