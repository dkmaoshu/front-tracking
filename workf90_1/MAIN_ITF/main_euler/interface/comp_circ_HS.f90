module compute_circulation_HS

use solution
! 'solution.f90' 

implicit none

public compute_circulation


contains


subroutine compute_circulation(circ)

implicit none

real(8), intent(out) :: circ

integer :: i, j
real(8) :: aaa
type(critical_cell), pointer :: tempg

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number)) 
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ndd(ndn))

call give_sth(0)
call give_cv
call give_gcl
call give_nd

circ=0.0d0

! Circulation along the top wall.
do i=nxll,nxll+nx-1
 aaa=x_velocity(uu(i,nyll+ny-1))
 circ=circ-aaa*h
end do

! Contribution along the horizontal centerline.
do i=nxll,nxll+nx-1
 select case(ggd_cell(i,0)%region)
  case('smth  ') 
   aaa=x_velocity(uu(i,0))
   circ=circ+aaa*h
  case('crit  ')
   call visit(ggd_cell(i,0)%ccpt(1)%address,tempg) 
   aaa=x_velocity(tempg%p_cell%l_state)
   circ=circ+aaa*h
  case default; call error_message
 end select
end do

! Contribution on the two ends
do j=nyll,nyll+ny-1
 aaa=y_velocity(uu(nxll,j))
 circ=circ-aaa*h
 aaa=y_velocity(uu(nxll+nx-1,j))
 circ=circ+aaa*h
end do    

call get_sth(1)
call get_cv
call get_nd
call get_gcl

deallocate(uu); deallocate(cvv)
deallocate(ggd_cell); deallocate(ndd)

end subroutine compute_circulation


end module compute_circulation_HS