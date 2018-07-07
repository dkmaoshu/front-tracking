module input_setup

use solution
! 'solution.f90'

implicit none

public input_and_setup
private uu, cvv, ggd_cell


contains


subroutine input_and_setup

implicit none

type(critical_cell), pointer :: current, temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: address, l_stk, r_stk
integer:: i, j, k, jj, nd_nmb, idx
type(curve_plug), pointer :: currentp, tempp
type(node_info) :: n_cell
type(cv_plug_info) :: plug
character*6 :: side_front, side_behind

open(4,file='d:\workf90_1\input\initt.dat') 
 read(4,'(f18.12)') current_time
close(4)
! Input the initial time.

open(1,file='d:\workf90_1\input\initg.dat')
 read(1,'(2f18.12)') x_width, y_width
 read(1,'(4i5)') nxll, nyll, nx, ny
close(1)
call set_grid
! Input the information of grid and set the grid.

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(cvn))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ndd(ndn))

do i=nxll-3,nxll+nx+2
 do j=nyll-3,nyll+ny+2
  uu(i,j)=error_data
 end do
end do

nxll_boundary=nxll; nyll_boundary=nyll
nx_boundary=nx; ny_boundary=ny

call set_sth
call set_gcl

open(1,file='d:\workf90_1\input\initu.dat')
read(1,'(4f20.12)') &
(((uu(i,j)%value(k),k=1,4),i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
close(1)
! Input the solution in smooth region.

open(5,file='d:\workf90_1\input\initgd.dat')
read(5,'(9a8)') &
((ggd_cell(i,j)%region, i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
read(5,'(15i5)') &
(((ggd_cell(i,j)%ccpt(jj)%address%cv_nmb, jj=1,4), &
i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
read(5,'(15i5)') &
(((ggd_cell(i,j)%ccpt(jj)%address%idx, jj=1,4), &
i=nxll,nxll+nx-1),j=nyll,nyll+ny-1)
read(5,'(8a8)') &
(((ggd_cell(i,j)%ccpt(jj)%side, jj=1,4), i=nxll,nxll+nx-1), &
j=nyll,nyll+ny-1) 
read(5,'(a2)') auxiliary_curve_type
read(5,'(a12)') order_of_curve_reset
close(5)
! Input critical cells' addresses on grid and pointers from
! grid cells to critical cells.

! call scan(ggd_cell)
! call scanning(uu)

open(2,file='d:\workf90_1\input\initd.dat')
do i=1,curves_number
 read(2,'(a8)') cvv(i)%status
 if(cvv(i)%status.eq.'awake'.or.cvv(i)%status.eq.'yawn  ') then 
  read(2,'(a8)') cvv(i)%cv_type
  read(2,'(2i5)') cvv(i)%begin_end, cvv(i)%end_end
  read(2,'(2i5)') cvv(i)%total, cvv(i)%wave
  allocate(cvv(i)%begin)
  nullify(cvv(i)%begin%previous); nullify(cvv(i)%begin%next)
  cvv(i)%begin%address%cv_nmb=i; cvv(i)%begin%address%idx=-1
    
  if(cvv(i)%status.eq.'awake ') then
   current=>cvv(i)%begin
   read(2,'(2i5)') address 
    
   do while(address%idx.gt.0) 
    
    read(2,'(i5,i5)') l_stk, r_stk
    
    read(2,'(a8)') g_cell%g_type
    read(2,'(2i5)') g_cell%x_idx, g_cell%y_idx
    read(2,'(2f18.12)') g_cell%dis_posi(1), g_cell%dis_posi(2)
    read(2,'(2i5)') g_cell%edge(1), g_cell%edge(2)
    read(2,'(4a8)') (g_cell%point(jj),jj=1,4)
    
    read(2,'(4f18.12)') (p_cell%l_state%value(k),k=1,4)
	read(2,'(4f18.12)') (p_cell%r_state%value(k),k=1,4)
	read(2,'(4f18.12)') (p_cell%or_state%value(k),k=1,4)
    read(2,'(i5)') p_cell%wv_nb
    
    allocate(temp)
     
    temp%address=address; temp%l_stk=l_stk; temp%r_stk=r_stk
    temp%g_cell=g_cell; temp%p_cell=p_cell
    current%next=>temp
    temp%previous=>current
    current=>temp
    read(2,'(2i5)') address
   end do
   
   if(cvv(i)%cv_type.eq.'circular') then
    current%next=>cvv(i)%begin%next
    cvv(i)%begin%next%previous=>current
   else
    nullify(cvv(i)%begin%next%previous)
   end if
  end if
   
  if(cvv(i)%status.eq.'yawn  ') read(2,'(2i5)') address
  allocate(cvv(i)%eend); nullify(cvv(i)%eend%next)
  cvv(i)%eend%address%cv_nmb=i; cvv(i)%eend%address%idx=-2
  if(associated(cvv(i)%begin%next)) then
   cvv(i)%eend%previous=>current
  else
   nullify(cvv(i)%eend%previous)
  end if

 end if
end do
close(2)
! Input discontinuity curves.

! call check_list_c(1,'down  ')

open(6,file='d:\workf90_1\input\initnd.dat')
do i=1,nodes_number
 read(6,'(a8)') ndd(i)%status
 if(ndd(i)%status.eq.'awake') then 
!  read(6,'(a8)') ndd(i)%n_type
  read(6,'(i5,i5)') n_cell%x_idx, n_cell%y_idx
  read(6,'(4f18.12)') (n_cell%or_state%value(k),k=1,4)
  ndd(i)%n_cell=n_cell

  read(6,'(2i5)') ndd(i)%cv_rg%total, ndd(i)%cv_rg%top_idx
  allocate(ndd(i)%cv_rg%begin)
  nullify(ndd(i)%cv_rg%begin%previous)
  nullify(ndd(i)%cv_rg%begin%next)
  ndd(i)%cv_rg%begin%address%idx=-1
  currentp=>ndd(i)%cv_rg%begin

  do while(.true.)
   read(6,'(2i5)') nd_nmb, idx 
   if(nd_nmb.le.0) exit
   read(6,'(2a8)') side_front, side_behind

   read(6,'(i5)') plug%cv_nmb
   read(6,'(a8)') plug%end_type
   read(6,'(i5)') plug%edge

   allocate(tempp)
   
   tempp%address%idx=idx; tempp%address%nd_nmb=nd_nmb
   tempp%plug=plug
   tempp%plug%side_front=side_front
   tempp%plug%side_behind=side_behind
   currentp%next=>tempp
   tempp%previous=>currentp
   currentp=>tempp
!   read(6,'(i5)') idx
  end do
  
  currentp%next=>ndd(i)%cv_rg%begin%next
  ndd(i)%cv_rg%begin%next%previous=>currentp

  allocate(ndd(i)%cv_rg%eend); nullify(ndd(i)%cv_rg%eend%next)
  ndd(i)%cv_rg%eend%address%idx=-2
  ndd(i)%cv_rg%eend%previous=>currentp

 end if
end do
close(6)
! Input discontinuity curves.

! call check_list_c(1,'up    ')
  
! call check_ring(1,'counter ')

! Input the boundary situation.
open(7,file='d:\workf90_1\input\boundary_situation.dat')
read(7,'(a12)') boundary_type
if(boundary_type.eq.'periodic    ') read(7,'(f18.12)') moved_time
close(7)

!call get_sth(0)
!call get_cv
!call get_gcl
!call get_nd

!deallocate(uu); deallocate(cvv)
!deallocate(ggd_cell); deallocate(ndd)

end subroutine input_and_setup


end module input_setup