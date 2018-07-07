module move_all
! This module moves all the region of solution by certain distance
! in certain direction.

use solution
! 'solution.f90'

implicit none

public  move_region
private move


contains


subroutine move_region(if_finishing)

implicit none

logical, intent(in) :: if_finishing
integer i0
real(8) :: time_difference

time_difference=current_time-moved_time
if(time_difference.ge.0.2d0.or.if_finishing) then
 i0=int(time_difference/h)

 call move(i0,i0)
 moved_time=moved_time+i0*h

!  print*, 'Moved cells=', i0
!  print*, 'Moved time =', moved_time
  
end if

end subroutine move_region


subroutine move(i0,j0)
! This subroutine moves the whole region of solution by `i' cells 
! in x-direction and `j' cells in y-direction.

implicit none
integer, intent(in) :: i0, j0

type(state), dimension(:,:), allocatable :: uu_t
type(grid_cell), dimension(:,:), allocatable :: ggd_t

integer :: i, j, l
type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
integer :: head_mark
logical :: head_switch

allocate(uu(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(cvv(cvn))
allocate(ggd_cell(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(ndd(ndn))
call give_sth(0)
call give_cv
call give_gcl
call give_nd

allocate(uu_t(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(ggd_t(nxll:nxll+nx-1,nyll:nyll+ny-1))

do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  uu_t(i,j)=0.0d0; call clean(ggd_t(i,j))
 end do
end do

do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  if(nxll.le.i+i0.and.nxll+nx-1.ge.i+i0) then
   if(nyll.le.j+j0.and.nyll+ny-1.ge.j+j0) then
    uu_t(i,j)=uu(i+i0,j+j0)
	ggd_t(i,j)=ggd_cell(i+i0,j+j0)
   end if
  end if    
 end do
end do

do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  uu(i,j)=uu_t(i,j); ggd_cell(i,j)=ggd_t(i,j)
 end do
end do
! The smooth region has been moved.

! call scanning(uu)

do l=1, cvn
 if(cvv(l)%status.eq.'awake ') then
  temp=>cvv(l)%begin%next
  head_mark=temp%address%idx; head_switch=.true.
  do while(associated(temp).and.head_switch)
   g_cell=temp%g_cell
   g_cell%x_idx=g_cell%x_idx-i0
   g_cell%y_idx=g_cell%y_idx-j0
   temp%g_cell=g_cell

   temp=>temp%next
   if(associated(temp)) then
    head_switch=(temp%address%idx.ne.head_mark)
   else
    exit
   end if
  end do
 end if
end do
! Discontinuity curves have been moved.

do l=1, ndn
 if(ndd(l)%status.eq.'awake ') then
  ndd(l)%n_cell%x_idx=ndd(l)%n_cell%x_idx-i0
  ndd(l)%n_cell%y_idx=ndd(l)%n_cell%y_idx-j0
 end if
end do

call get_sth(0)
call get_cv
call get_gcl
call get_nd

deallocate(uu); deallocate(cvv)
deallocate(ggd_cell); deallocate(ndd)
deallocate(uu_t); deallocate(ggd_t)

end subroutine move


end module move_all