module grid_map

use grid
! 'grid.f90'

use adss_info_cell
! 'adss_cell.f90'

use discontinuity_curves
! 'discontinuity_curve.f90'

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use solu_in_smth
! 'solution_smooth.f90'

implicit none

type corner_pointer
 type(adss_info) :: address
 character*6 :: side
end type corner_pointer

type grid_cell
 character*6 :: region
 type(corner_pointer), dimension(4) :: ccpt
end type grid_cell
! The component region tells whether the cell is occupied by 
! critical cells, a node cell, or in smooth region. The component
! ccpt is pointers at the four corners pointing to the critical
! cell that occupies the cell.

type(grid_cell), dimension(:,:), allocatable :: gd_cell
type(grid_cell), dimension(:,:), allocatable :: ggd_cell

interface assignment(=)
 module procedure assign_p
end interface

interface scanning
 module procedure scan_region
end interface

interface clean
 module procedure clean_gd, clean_ccpt
end interface

interface remove
 module procedure remove_c, remove_ac
end interface

interface operator(.eq.)
 module procedure compare_truee
end interface

interface operator(.ne.)
 module procedure compare_falsee
end interface

private gd_cell, scan_region, clean_gd, remove_c, compare_truee, compare_falsee
public  ggd_cell, set_gcl, give_gcl, get_gcl, scanning,  &
        clean_up_grid_map, clean, remove, update_map


contains


subroutine assign_p(ccpt1,ccpt2)

implicit none
type(corner_pointer), intent(out) :: ccpt1
type(corner_pointer), intent(in) :: ccpt2

ccpt1%side=ccpt2%side
ccpt1%address=ccpt2%address

end subroutine assign_p


subroutine set_gcl

implicit none

allocate(gd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))

end subroutine set_gcl


subroutine give_gcl

implicit none
integer :: i, j

ggd_cell%region=gd_cell%region
do i=nxll-3,nxll+nx+2
 do j=nyll-3, nyll+ny+2
  ggd_cell(i,j)%ccpt=gd_cell(i,j)%ccpt
 end do
end do

end subroutine give_gcl


subroutine get_gcl

implicit none
integer :: i, j

gd_cell%region=ggd_cell%region
do i=nxll-3,nxll+nx+2
 do j=nyll-3,nyll+ny+2
  gd_cell(i,j)%ccpt=ggd_cell(i,j)%ccpt
 end do
end do

end subroutine get_gcl


subroutine clean_up_grid_map(gd_cell)
! Clean-up grid map in all grid cells.

implicit none

type(grid_cell), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2) :: gd_cell
integer :: i, j

do i=nxll-3,nxll+nx+2
 do j=nyll-3,nyll+ny+2
  call clean(gd_cell(i,j))
 end do
end do

end subroutine clean_up_grid_map


subroutine clean_gd(gdx_cell)

implicit none

type(grid_cell), intent(out) ::	gdx_cell

gdx_cell%region='smth  '
gdx_cell%ccpt%address=adss_info(0,0)
gdx_cell%ccpt%side=' '

end subroutine clean_gd


subroutine scan_region(gd_cell)

implicit none

type(grid_cell), dimension(:,:) :: gd_cell

type(grid_cell), dimension(-3:3,-3:3) :: reg
integer :: i0, j0, i, j

write(*,*) ' Please input I0 and J0.'
read(*,'(2i5)') i0, j0

reg=gd_cell(i0-nxll-2:i0-nxll+4,j0-nyll-2:j0-nyll+4)
do j=3,-3,-1
 write(*,'(7a10)') (reg(i,j)%region, i=-3,3)
end do

end subroutine scan_region

subroutine find_neighbor_cell(i_c,j_c,across_edge,i_n,j_n)
! This subroutine finds the neighbor cell of a grid cell across
! a given edge.

implicit none
integer, intent(in) :: i_c, j_c
! The indexes of the given grid cell
integer, intent(in) :: across_edge
! The given across edge
integer, intent(out) :: i_n,j_n
! The indexes of the neighbor grid cell.

i_n=i_c; j_n=j_c
select case(across_edge)
 case(1); j_n=j_n-1
 case(2); i_n=i_n+1
 case(3); j_n=j_n+1
 case(4); i_n=i_n-1
end select

end subroutine find_neighbor_cell


subroutine remove_c(temp,position,remove_side,difference)
! This subroutine removes a critical cell from a discontinuity
! curve list. 
! The difference between 'remove' here and 'delete' defined in 
! module 'discontinuity_curves' is that 'delete' only deletes
! the critical cell, but 'remove' deletes the critical cell, 
! updates the space neighboring relation and grid map, and picks
! the conservation-differnce caused by the remove.

implicit none
type(critical_cell), pointer :: temp
character*6, intent(in) :: position
! For the meaning of the above two variables, see 'delete' in
! module 'discontinuity_curves'.
character*6, intent(in) :: remove_side
! Tell from which side the critical cell is to be removed.
type(state), intent(out) :: difference
! The conservation difference caused by the mremove.

type(critical_cell), pointer :: temp_l, temp_r
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: adss, adss_l, adss_r
character*6 :: from_side
integer :: i, i0, j0

nullify(temp_l); nullify(temp_r)
difference=error_data
g_cell=temp%g_cell; p_cell=temp%p_cell
i0=g_cell%x_idx; j0=g_cell%y_idx

! The following updates the space neighboring relation and the 
! grid map.
adss=temp%address; adss_l=temp%l_stk; adss_r=temp%r_stk
if(adss_l%cv_nmb.le.0.and.adss_r%cv_nmb.le.0) then
 call clean(ggd_cell(i0,j0))
end if

if(adss_l%cv_nmb.le.0.and.adss_r%cv_nmb.gt.0) then
 call visit(adss_r,temp_r)
 if(temp_r%l_stk.eq.adss) then
  from_side='left  '; temp_r%l_stk=adss_info(0,0)
 else
  if(temp_r%r_stk.eq.adss) then
   from_side='right '; temp_r%r_stk=adss_info(0,0)
  else
   print*, 'Something wrong'; pause
  end if
 end if
 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%address.eq.adss) then
   ggd_cell(i0,j0)%ccpt(i)%address=adss_r
   ggd_cell(i0,j0)%ccpt(i)%side=from_side
  end if
 end do
end if  

if(adss_l%cv_nmb.gt.0.and.adss_r%cv_nmb.le.0) then
 call visit(adss_l,temp_l)
 if(temp_l%l_stk.eq.adss) then
  from_side='left  '; temp_l%l_stk=adss_info(0,0)
 else
  if(temp_l%r_stk.eq.adss) then
   from_side='right '; temp_l%r_stk=adss_info(0,0)
  else
   print*, 'Something wrong'; pause
  end if
 end if
 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%address.eq.adss) then
   ggd_cell(i0,j0)%ccpt(i)%address=adss_l
   ggd_cell(i0,j0)%ccpt(i)%side=from_side
  end if
 end do
end if

if(adss_l%cv_nmb.gt.0.and.adss_r%cv_nmb.gt.0) then
 call visit(adss_l,temp_l); call visit(adss_r,temp_r)
 if(temp_l%l_stk.eq.adss) then
  temp_l%l_stk=adss_r
 else
  if(temp_l%r_stk.eq.adss) then
   temp_l%r_stk=adss_r
  else
   print*, 'Something wrong'; pause
  end if
 end if
 if(temp_r%l_stk.eq.adss) then
  temp_r%l_stk=adss_l
 else
  if(temp_r%r_stk.eq.adss) then
   temp_r%r_stk=adss_l
  else
   print*, 'Something wrong'; pause
  end if
 end if
end if

! The following picks the conservation-difference
select case(remove_side)
 case('left  ')
  difference=p_cell%or_state-p_cell%r_state
  uu(i0,j0)=p_cell%r_state
 case('right ')
  difference=p_cell%or_state-p_cell%l_state
  uu(i0,j0)=p_cell%l_state
 case('middle')
  difference=0.d0
end select

! Delete the critical cell from the discontinuity curve list.
call deletee(temp,position)

end subroutine remove_c

subroutine remove_ac(temp,position,remove_side,difference)
! This subroutine removes a critical cell from a discontinuity
! curve list. 
! The difference between 'remove' here and 'delete' defined in 
! module 'discontinuity_curves' is that 'delete' only deletes
! the critical cell, but 'remove' deletes the critical cell, 
! updates the space neighboring relation and grid map, and picks
! the conservation-differnce caused by the remove.

implicit none
type(auxi_crit_cell), pointer :: temp
character*6, intent(in) :: position
! For the meaning of the above two variables, see 'delete' in
! module 'discontinuity_curves'.
character*6, intent(in) :: remove_side
! Tell from which side the critical cell is to be removed.
type(state), intent(out) :: difference
! The conservation difference caused by the mremove.

type(auxi_crit_cell), pointer :: temp_l, temp_r
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: adss, adss_l, adss_r
character*6 :: from_side
integer :: i, i0, j0

nullify(temp_l); nullify(temp_r)
difference=error_data
g_cell=temp%g_cell; p_cell=temp%p_cell
i0=g_cell%x_idx; j0=g_cell%y_idx

! The following updates the space neighboring relation and the 
! grid map.
adss=temp%address; adss_l=temp%l_stk; adss_r=temp%r_stk
if(adss_l%cv_nmb.le.0.and.adss_r%cv_nmb.le.0) then
 call clean(ggd_cell(i0,j0))
end if

if(adss_l%cv_nmb.le.0.and.adss_r%cv_nmb.gt.0) then
 call visit(adss_r,temp_r)
 if(temp_r%l_stk.eq.adss) then
  from_side='left  '; temp_r%l_stk=adss_info(0,0)
 else
  if(temp_r%r_stk.eq.adss) then
   from_side='right '; temp_r%r_stk=adss_info(0,0)
  else
   print*, 'Something wrong'; pause
  end if
 end if
 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%address.eq.adss) then
   ggd_cell(i0,j0)%ccpt(i)%address=adss_r
   ggd_cell(i0,j0)%ccpt(i)%side=from_side
  end if
 end do
end if  

if(adss_l%cv_nmb.gt.0.and.adss_r%cv_nmb.le.0) then
 call visit(adss_l,temp_l)
 if(temp_l%l_stk.eq.adss) then
  from_side='left  '; temp_l%l_stk=adss_info(0,0)
 else
  if(temp_l%r_stk.eq.adss) then
   from_side='right '; temp_l%r_stk=adss_info(0,0)
  else
   print*, 'Something wrong'; pause
  end if
 end if
 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%address.eq.adss) then
   ggd_cell(i0,j0)%ccpt(i)%address=adss_l
   ggd_cell(i0,j0)%ccpt(i)%side=from_side
  end if
 end do
end if

if(adss_l%cv_nmb.gt.0.and.adss_r%cv_nmb.gt.0) then
 call visit(adss_l,temp_l); call visit(adss_r,temp_r)
 if(temp_l%l_stk.eq.adss) then
  temp_l%l_stk=adss_r
 else
  if(temp_l%r_stk.eq.adss) then
   temp_l%r_stk=adss_r
  else
   print*, 'Something wrong'; pause
  end if
 end if
 if(temp_r%l_stk.eq.adss) then
  temp_r%l_stk=adss_l
 else
  if(temp_r%r_stk.eq.adss) then
   temp_r%r_stk=adss_l
  else
   print*, 'Something wrong'; pause
  end if
 end if
end if

! The following picks the conservation-difference
select case(remove_side)
 case('left  ')
  difference=p_cell%or_state-p_cell%r_state
  uu(i0,j0)=p_cell%r_state
 case('right ')
  difference=p_cell%or_state-p_cell%l_state
  uu(i0,j0)=p_cell%l_state
 case('middle')
  difference=0.d0
end select

! Delete the critical cell from the discontinuity curve list.
call deletee(temp,position)

end subroutine remove_ac


subroutine update_map(t_new,side)
! This subroutine updates the grid-map in a grid cell when a
! newly produced critical cell is inserted with no other critical
! cell stacked in the same cell on the 'side' side.

implicit none
type(critical_cell), pointer :: t_new
character*6, intent(in) :: side

type(geo_info) :: g_cell
integer :: i0, j0, i

g_cell=t_new%g_cell; i0=g_cell%x_idx; j0=g_cell%y_idx
ggd_cell(i0,j0)%region='crit  '
do i=1, 4
 if(g_cell%point(i).eq.side) then
  ggd_cell(i0,j0)%ccpt(i)%address=t_new%address
  ggd_cell(i0,j0)%ccpt(i)%side=side
 end if
end do

end subroutine update_map


function compare_truee(ccpt_1,ccpt_2) result(c) 

implicit none
type(corner_pointer), intent(in) :: ccpt_1, ccpt_2
logical :: c

c=.false.
if(ccpt_1%address.eq.ccpt_2%address) then
 if(ccpt_1%side.eq.ccpt_2%side) c=.true.
end if

end function compare_truee


function compare_falsee(ccpt_1,ccpt_2) result(c) 

implicit none
type(corner_pointer), intent(in) :: ccpt_1, ccpt_2
logical :: c

c=.false.
if(ccpt_1%address.ne.ccpt_2%address) then
 c=.true.
else
 if(ccpt_1%side.eq.ccpt_2%side) c=.true.
end if

end function compare_falsee


subroutine clean_ccpt(ccpt)

implicit none
type(corner_pointer), intent(out) :: ccpt

call clean_up_address(ccpt%address); ccpt%side='      '

end  subroutine clean_ccpt


end module grid_map