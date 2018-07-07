module out_intermediate
! This module describes the operation of output of intermediate
! results.

use solution
! 'solution.f90'

use grid_map
! 'grid_sth.f90'

implicit none


private uu, cvv, ggd_cell
public output_m


contains


subroutine output_m

implicit none

type(critical_cell), pointer :: current, temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: address, l_stk, r_stk
integer:: head_mark, i, j, k, jj, nd_nmb, idx
logical :: head_switch
type(curve_plug), pointer :: tempp, currentp
type(cv_plug_info) :: plug
type(node_info) :: n_cell
character*6 :: side_front, side_behind
type(state) :: fff

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ndd(nodes_number))
call give_sth(0)
call give_cv
call give_gcl
call give_nd

open(4,file='d:\workf90_1\output\initt.dat')
write(4,'(f18.12)') current_time
write(4,'(i5)') total_step
close(4)
! Output the final time.

open(3,file='d:\workf90_1\output\initg.dat')
write(3,'(2f18.12)') x_width, y_width
write(3,'(4i5)') nxll, nyll, nx, ny
close(3)
! Output the information of grid.

!call scan_sth(uu)

open(1,file='d:\workf90_1\output\initu.dat')
write(1,'(4f20.12)')  &
(((uu(i,j)%value(k),k=1,4),i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
write(1,'(4f20.12)')  &
((uu(i,j)%gamma,i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)

close(1) 
! Output the solution in smooth region.

open(5,file='d:\workf90_1\output\initgd.dat')
write(5,'(9a8)') &
((ggd_cell(i,j)%region, i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
write(5,'(15i5)') &
(((ggd_cell(i,j)%ccpt(jj)%address%cv_nmb, jj=1,4), &
i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
write(5,'(15i5)') &
(((ggd_cell(i,j)%ccpt(jj)%address%idx, jj=1,4), &
i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
write(5,'(8a8)') &
(((ggd_cell(i,j)%ccpt(jj)%side, jj=1,4), i=nxll,nxll+nx-1), &
j=nyll,nyll+ny-1)
write(5,'(a2)') auxiliary_curve_type
write(5,'(a12)') order_of_curve_reset
close(5)
! Output critical cells' addresses on grid and pointers from
! grid cells to critical cells.

open(2,file='d:\workf90_1\output\initd.dat')
do i=1,curves_number
 write(2,'(a8)') cvv(i)%status
 if(cvv(i)%status.eq.'awake'.or.cvv(i)%status.eq.'yawn  ') then
  write(2,'(a8)') cvv(i)%cv_type
  write(2,'(2i5)') cvv(i)%begin_end, cvv(i)%end_end
  write(2,'(2i5)') cvv(i)%total, cvv(i)%wave
  current=>cvv(i)%begin; head_switch=.true.
  if(associated(current%next)) head_mark=current%next%address%idx
  do while(associated(current%next).and.head_switch)

   temp=>current%next
   g_cell=temp%g_cell; p_cell=temp%p_cell
   address=temp%address; l_stk=temp%l_stk; r_stk=temp%r_stk

   fff=p_cell%l_state
   fff=p_cell%r_state
   fff=p_cell%or_state

   write(2,'(i5,i5)') address
   write(2,'(i5,i5)') l_stk, r_stk

   write(2,'(a8)') g_cell%g_type
   write(2,'(2i5)') g_cell%x_idx, g_cell%y_idx
   write(2,'(2f18.12)') g_cell%dis_posi(1), g_cell%dis_posi(2)
   write(2,'(2i5)') g_cell%edge(1), g_cell%edge(2)
   write(2,'(4a8)') (g_cell%point(jj),jj=1,4)
   
   write(2,'(4f18.12)') (p_cell%l_state%value(k),k=1,4)
   write(2,'(f20.12)') p_cell%l_state%gamma
   write(2,'(4f18.12)') (p_cell%r_state%value(k),k=1,4)
   write(2,'(f20.12)') p_cell%r_state%gamma
   write(2,'(4f18.12)') (p_cell%or_state%value(k),k=1,4)
   write(2,'(f20.12)') p_cell%or_state%gamma
   write(2,'(i5)') p_cell%wv_nb

   current=>temp
   if(associated(temp%next))  &
   head_switch=(temp%next%address%idx.ne.head_mark)
  end do
  write(2,'(i5,i5)') i, -1
 end if
end do
close(2)
! Output discontinuity curves.

open(6,file='d:\workf90_1\output\initnd.dat')
do i=1,nodes_number
 write(6,'(a8)') ndd(i)%status
 if(ndd(i)%status.eq.'awake') then
  n_cell=ndd(i)%n_cell
!  write(6,'(a8)') ndd(i)%n_type
  write(6,'(i5,i5)') n_cell%x_idx, n_cell%y_idx
  write(6,'(4f18.12)') (n_cell%or_state%value(k),k=1,4)

  write(6,'(2i5)') ndd(i)%cv_rg%total, ndd(i)%cv_rg%top_idx
  currentp=>ndd(i)%cv_rg%begin; head_switch=.true.
  head_mark=currentp%next%address%idx
  do while(head_switch)
   tempp=>currentp%next
   nd_nmb=tempp%address%nd_nmb; idx=tempp%address%idx
   side_front=tempp%plug%side_front 
   side_behind=tempp%plug%side_behind
   plug=tempp%plug

   write(6,'(2i5)') nd_nmb, idx
   write(6,'(2a8)') side_front, side_behind
   write(6,'(i5)') plug%cv_nmb
   write(6,'(a8)') plug%end_type
   write(6,'(i5)') plug%edge

   currentp=>tempp
   head_switch=(tempp%next%address%idx.ne.head_mark)
  end do
  write(6,'(2i5)') 0, -1
 end if
end do
close(6)
! Output node cells information.

! Output the boundary situation.
open(7,file='d:\workf90_1\output\boundary_situation.dat')
write(7,'(a12)') boundary_type
if(boundary_type.eq.'periodic    ') write(7,'(f18.12)') moved_time
close(7)


deallocate(uu); deallocate(cvv) 
deallocate(ggd_cell); deallocate(ndd)

end subroutine output_m

end module out_intermediate