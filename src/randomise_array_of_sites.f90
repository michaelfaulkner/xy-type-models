subroutine randomise_array_of_sites
use variables
implicit none
integer :: i, j, store_site

do i = 1, no_of_sites
    j = 1 + floor(no_of_sites * rand())
    store_site = array_of_sites(i)
    array_of_sites(i) = array_of_sites(j)
    array_of_sites(j) = store_site
end do

return
end subroutine randomise_array_of_sites
