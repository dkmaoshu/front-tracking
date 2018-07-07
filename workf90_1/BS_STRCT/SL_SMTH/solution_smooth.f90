module solu_in_smth

use physical_state
! 'burgers.f90

use grid
! 'grid.f90'

use show_smth
! 'sh_sths.f90'

type(state), dimension(:,:), allocatable :: u, u1
type(state), dimension(:,:), allocatable :: uu

interface scanning
 module procedure scan_sth
end interface

public  uu, set_sth, give_sth, get_sth, scanning
private u, u1, scan_sth


contains


subroutine set_sth

implicit none
integer :: i, j

allocate(u(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(u1(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
do i=nxll-3, nxll+nx+2
 do j=nyll-3, nyll+ny+2
  u1(i,j)=error_data
 end do
end do

end subroutine set_sth


subroutine give_sth(level)

implicit none
integer :: level

if(level.eq.0) then
 uu=u
else
 uu=u1
end if

end subroutine give_sth


subroutine get_sth(level)

implicit none
integer :: level

if(level.eq.0) then
 u=uu
else
 u1=uu
end if

end subroutine get_sth


subroutine scan_sth(u)
implicit none
type(state), dimension(:,:) :: u

type(state), dimension(-3:3,-3:3) :: u_sh
integer :: i0, j0, i, j

write(*,*) ' Please input I0 and J0.'
read(*,'(2i5)') i0, j0

do i=-3,3
 do j=-3,3
  u_sh(i,j)=u(i+i0-nxll+1,j+j0-nyll+1)
 end do
end do
call show(u_sh)

end subroutine scan_sth


end module solu_in_smth