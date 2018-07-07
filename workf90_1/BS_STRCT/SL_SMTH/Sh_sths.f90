module show_smth

use physical_state
implicit none


contains


subroutine show(u_sh)

implicit none

type(state), dimension(-3:3,-3:3) :: u_sh
integer :: i, j

do j=3,-3,-1
 write(*,'(7f10.4)') (u_sh(i,j)%value,i=-3,3)
end do

end subroutine show

end module show_smth