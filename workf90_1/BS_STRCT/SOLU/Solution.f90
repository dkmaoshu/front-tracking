module solution

use grid_map
! 'grid_map.f90'

use solu_in_smth
! 'solution_smooth.f90'

use discontinuity_curves
! 'discontinuity_curve.f90'

use node_cells
! 'nd_cls.f90'

implicit none

type crit_info_c
 type(adss_info) :: address
 type(geo_info) :: g_cell
 type(phy_info) :: p_cell
 type(state), dimension(0:2) :: l_state, r_state, or_state
 character*6 :: side
end type crit_info_c

type(crit_info_c) :: crit

interface smth_dpr
 module procedure smth_dpr_c
end interface

public  uu, cvv, ggd_cell, smth_dpr, find_other_end, &
        evaluation_in_cell, change_node_number
private pick_go_ord, pick_go_crt, smth_dpr_c, transf, crit, &
        pick_smth


contains


subroutine smth_dpr_c(i,j,temp,uxy,xy_dir,side,start,if_expand,hg_ord,range_slide)
! This subroutine prepares the data for operations in smooth 
! regions. At grid cells that belong to the smooth region, it 
! maintains the original data of the solution, and at grid cells 
! that are out of the smooth region, it supplies the extrapolated	
! data from the data of the solution within the region. The 
! variable 'hg_ord' is the highest order of extrapolation +1 to 
! be used.

implicit none

integer, intent(in) :: i, j
! The indexes of the starting grid cell.
type(critical_cell), pointer :: temp
! Containing the information of the starting critical cell.
type(state), dimension(-3:3) :: uxy
! The data sector.
character*3, intent(in) :: xy_dir
! The direction of the preparing.
character*6, intent(in) :: side
character*8, intent(in) :: start
! Indicating whether starting with a regular or critical cell.
integer, intent(in) :: hg_ord
! The highest order of extrapolation+1.
logical, intent(in) :: if_expand
! Determine the range of the to-be expanded data, the possible values are 'forward ', 'backward'
! and 'full    '.
character*8, intent(in) :: range_slide 

integer :: kkk
! real(8), dimension(2) :: pt

call transf(temp,side)
kkk=i
kkk=j
call pick_smth(i,j,uxy,xy_dir,start)
if(if_expand) call data_fill(uxy,hg_ord,range_slide)

! call find_corner(pt,1)

end subroutine smth_dpr_c


subroutine pick_smth(i,j,uxy,xy_dir,start)
! This subroutine picks up the smooth data starting either from
! an ordinary cell or a side of a critical cell.

implicit none

integer, intent(in) :: i, j
type(state), dimension(-3:3) :: uxy
character*3, intent(in) :: xy_dir
character*8, intent(in) :: start
! Check on subroutine 'sthm_dpr' for the meaning of the above 
! variables.

type(crit_info_c) :: cr
integer :: i0, j0
character*8 :: stoped
integer :: k, kl

do kl=-3,3
 uxy(kl)=error_data
end do
! Initialization.

stoped=start
if(start.eq.'ordinary') then
 i0=i; j0=j
else
 cr=crit
end if
do k=0,3
 if(stoped.eq.'ordinary') then
  call clean_up_g_cell(cr%g_cell); call clean_up_p_cell(cr%p_cell)
  cr%address=adss_info(0,0); cr%side='     '
  if(xy_dir.eq.'xx') then
   call pick_go_ord(i0,j0,cr,uxy(k),1,stoped)
  else
   call pick_go_ord(i0,j0,cr,uxy(k),2,stoped)
  end if
 else 
  if(stoped.eq.'critical') then
   if(xy_dir.eq.'xx') then
    call pick_go_crt(i0,j0,cr,uxy(k),1,stoped)
   else
    call pick_go_crt(i0,j0,cr,uxy(k),2,stoped)
   end if
  else
   exit
  end if
 end if
end do
! Pick-up in positive direction.

stoped=start
if(start.eq.'ordinary') then
 i0=i; j0=j
else
 cr=crit
end if
do k=0,-3,-1
 if(stoped.eq.'ordinary') then
  call clean_up_g_cell(cr%g_cell); call clean_up_p_cell(cr%p_cell)
  cr%address=adss_info(0,0); cr%side='     '
  if(xy_dir.eq.'xx') then
   call pick_go_ord(i0,j0,cr,uxy(k),3,stoped)
  else
   call pick_go_ord(i0,j0,cr,uxy(k),4,stoped)
  end if
 else
  if(stoped.eq.'critical') then
   if(xy_dir.eq.'xx') then
    call pick_go_crt(i0,j0,cr,uxy(k),3,stoped)
   else
    call pick_go_crt(i0,j0,cr,uxy(k),4,stoped)
   end if
  else
   exit
  end if
 end if
end do
! Pick-up in negative direction.

end subroutine pick_smth


subroutine pick_go_ord(i0,j0,cr,value,direction,stoped)
! This subroutine performces the pick-and-go operation starting 
! with an ordinary cell to collect smooth data. 

implicit none

integer, intent(inout) :: i0, j0
! Indexes of the ordinary cell started with and also of the 
! possible stop-at ordinary cell.
type(crit_info_c), intent(out) :: cr
! Possible stop-at critical cell and its side.
type(state), intent(out) :: value
! Picked smooth data in the starting ordinary cell.
integer, intent(in) :: direction
! In which directon the GO is performed, the possible values
! are 1, 2, 3 and 4 indicating the x-positive, y-positive, 
! x-negative and y-negative directions, respectively.
character*8, intent(out) :: stoped
! Indicating the kind of the cell stoped at, the possible values
! are 'ordinary' and 'critical'.

type(critical_cell), pointer :: temp
integer :: chosen_corner

value=uu(i0,j0)
! Pick up the data in the starting ordinary cell.

select case(direction)
 case(1); i0=i0+1
 case(2); j0=j0+1
 case(3); i0=i0-1
 case(4); j0=j0-1
end select
! Make a move in the direction.

if(i0.lt.nxll.or.i0.gt.nxll+nx-1) stoped='blind'
if(j0.lt.nyll.or.j0.gt.nyll+ny-1) stoped='blind'
if(stoped.ne.'blind') then
 if(ggd_cell(i0,j0)%region.eq.'smth') then
  stoped='ordinary'
 else
  if(ggd_cell(i0,j0)%region.eq.'crit') then
   stoped='critical'
   chosen_corner=cycle_4(direction)
   cr%address=ggd_cell(i0,j0)%ccpt(chosen_corner)%address
   call visit(cr%address,temp)
   cr%g_cell=temp%g_cell; cr%p_cell=temp%p_cell
   cr%side=ggd_cell(i0,j0)%ccpt(chosen_corner)%side
  else
   if(ggd_cell(i0,j0)%region.eq.'node') then
    stoped='blind'
   else
    print*, 'There must be something wrong in pick_go_ord.'
    stop
   end if           
  end if
 end if
end if   
! Identify the stoped-at cell.

end subroutine pick_go_ord


subroutine pick_go_crt(i0,j0,cr,value,direction,stoped)
! This subroutine performs the pick-and-go operation starting 
! with a critical cell to collect smooth data. 

implicit none

integer, intent(out) :: i0, j0
! Indexes of the possible stop-at ordinary cell.
type(crit_info_c), intent(inout) :: cr
! The critical cell started with.
type(state), intent(out) :: value
! Picked smooth data in the starting critical cell.
integer, intent(in) :: direction
! In which directon the GO is performed, the possible values
! are 1, 2, 3 and 4 indicating the x-positive, y-positive, 
! x-negative and y-negative directions, respectively.
character*8, intent(out) :: stoped
! Indicating the kind of the cell stoped at, the possible values
! are 'ordinary', 'critical' and 'blind', the last one means 
! that no further go is allowed.

type(critical_cell), pointer :: temp, temp_c, temp_n
integer, dimension(2) :: corners
integer :: chosen_corner, k
character*3 :: dir
type(geo_info) :: gg_cl
type(phy_info) :: pp_cl

type(corner_pointer) :: right_ccpt, also_right_ccpt
! Critical cells may stacked in the same grid cell. When this
! happens, 'corner_pointer's at grid points between the two
! stacked critical cells can be valued in two different ways as
! viewed from the two critical cells. The 'right_ccpt' and the
! 'also_right_ccpt' are these two values, and their used to 
! determine whether we should stop or go further to pick smooth 
! data of the solution.

i0=cr%g_cell%x_idx; j0=cr%g_cell%y_idx
if(cr%side.eq.'left') then
 value=cr%p_cell%l_state
else
 if(cr%side.eq.'right') then
  value=cr%p_cell%r_state
 else
  print*, 'There must be something wrong in pick_go_crit'
  stop
 end if
end if
! Pick up the data in the starting ordinary cell.

! Evaluate the 'right_ccpt' and 'also_right_ccpt'.
call visit(cr%address,temp_c)
right_ccpt%address=temp_c%address; right_ccpt%side=cr%side

call clean(also_right_ccpt); nullify(temp_n)
select case(cr%side)
 case('left  ')
  if(temp_c%l_stk%cv_nmb.gt.0) call visit(temp_c%l_stk,temp_n)
 case('right ')
  if(temp_c%r_stk%cv_nmb.gt.0) call visit(temp_c%r_stk,temp_n)
end select

also_right_ccpt%address%cv_nmb=error_data
also_right_ccpt%address%idx=error_data


if(associated(temp_n)) then
 also_right_ccpt%address=temp_n%address
 if(temp_n%l_stk.eq.temp_c%address) then
  also_right_ccpt%side='left  '
 else
  if(temp_n%r_stk.eq.temp_c%address) then
   also_right_ccpt%side='right '
  else
   call error_message
  end if
 end if
end if

stoped='        '
do k=1,2; corners(k)=cycle_4(direction+k); end do

if(ggd_cell(i0,j0)%ccpt(corners(1)).eq.right_ccpt.or. &
ggd_cell(i0,j0)%ccpt(corners(1)).eq.also_right_ccpt) then
 chosen_corner=corners(1)
else
 if(ggd_cell(i0,j0)%ccpt(corners(2)).eq.right_ccpt.or. &
 ggd_cell(i0,j0)%ccpt(corners(2)).eq.also_right_ccpt) then
  chosen_corner=corners(2)
 else
  stoped='blind'
 end if
end if

dir='   '
if(stoped.ne.'blind') then  
 select case(direction)
  case(1); i0=i0+1; dir='xx '
  case(2); j0=j0+1; dir='yy '
  case(3); i0=i0-1; dir='xx '
  case(4); j0=j0-1; dir='yy '
 end select
 chosen_corner=point_shift(chosen_corner,dir)
end if
! Make a move in the direction.

if(i0.lt.nxll.or.i0.gt.nxll+nx-1) stoped='blind'
if(j0.lt.nyll.or.j0.gt.nyll+ny-1) stoped='blind'
if(stoped.ne.'blind') then
 if(ggd_cell(i0,j0)%region.eq.'smth') then
  stoped='ordinary'
 else
  if(ggd_cell(i0,j0)%region.eq.'crit') then
   stoped='critical'
   cr%address=ggd_cell(i0,j0)%ccpt(chosen_corner)%address
   call visit(cr%address,temp)
   cr%g_cell=temp%g_cell; cr%p_cell=temp%p_cell
   cr%side=ggd_cell(i0,j0)%ccpt(chosen_corner)%side
  else
   if(ggd_cell(i0,j0)%region.eq.'node') then
    stoped='blind'
   else
    print*, 'There must be something wrong in pick_go_crt.'
    stop
   end if
  end if
 end if
end if
! Identify the stoped-at cell.

! The following patch is to deal with the case that two critical
! cells on different curves are stacked in the same grid cell
! with no corner point of the grid cell between them
if(stoped.eq.'blind   ') then
 
 chosen_corner=error_index
 do k=1,2
  if(cr%g_cell%point(corners(k)).eq.cr%side) then
   chosen_corner=corners(k); exit
  end if
 end do
 
 if(chosen_corner.gt.error_index) then
 
  select case(direction)
   case(1); i0=i0+1; dir='xx '
   case(2); j0=j0+1; dir='yy '
   case(3); i0=i0-1; dir='xx '
   case(4); j0=j0-1; dir='yy '
  end select

  call visit(cr%address,temp)
  if(associated(temp%previous)) then
   temp_n=>temp%previous
   gg_cl=temp_n%g_cell; pp_cl=temp_n%p_cell
   if(i0.eq.gg_cl%x_idx.and.j0.eq.gg_cl%y_idx) then
    cr%g_cell=gg_cl; cr%p_cell=pp_cl
    stoped='critical'
   end if
  end if
  if(stoped.eq.'blind   ') then
   if(associated(temp%next)) then
    temp_n=>temp%next
    gg_cl=temp_n%g_cell; pp_cl=temp_n%p_cell
    if(i0.eq.gg_cl%x_idx.and.j0.eq.gg_cl%y_idx) then
     cr%g_cell=gg_cl; cr%p_cell=pp_cl
     stoped='critical'
    end if
   end if
  end if
    
 end if  
endif   

end subroutine pick_go_crt


subroutine transf(temp,side)

implicit none

type(critical_cell), pointer :: temp
character*6, intent(in) :: side
type(critical_cell) :: fff

integer :: i

call clean_up_address(crit%address)
call clean_up_g_cell(crit%g_cell); call clean_up_p_cell(crit%p_cell)
crit%side='      '
do i=0,2
 crit%l_state(i)=error_data; crit%r_state(i)=error_data
 crit%or_state(i)=error_data
end do

if(associated(temp)) then
 crit%address=temp%address
 fff%address=crit%address
 crit%g_cell=temp%g_cell; crit%p_cell=temp%p_cell
 crit%side=side
end if

end subroutine transf


subroutine find_other_end(this_tempp,other_tempp)
! Given a curve-plug, this subroutine finds the other end of the
! curve.

implicit none
type(curve_plug), pointer :: this_tempp, other_tempp

type(cv_plug_info) :: plug
integer :: cv_nmb, nd_nmb
integer :: head_mark
logical :: head_switch

! Find the node cell the other end plugs in.
plug=this_tempp%plug; cv_nmb=plug%cv_nmb
select case(plug%end_type)
 case('begin '); nd_nmb=cvv(cv_nmb)%end_end
 case('end   '); nd_nmb=cvv(cv_nmb)%begin_end
 case default; print*, 'something is wrong!!!'; pause
end select

! Find the curve-plug of the other end.
other_tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=other_tempp%address%idx; head_switch=.true.
do while(head_switch)
 plug=other_tempp%plug
 if(plug%cv_nmb.eq.cv_nmb) exit
 other_tempp=>other_tempp%next
 head_switch=(other_tempp%address%idx.ne.head_mark)
end do
if(.not.associated(other_tempp)) then
 print*, 'Something is wrong!!!'; pause
end if

end subroutine find_other_end


subroutine change_node_number(old_num,new_num)
! This subroutine changes the node number 'old_num' to 'new_num'.

implicit none
integer, intent(in) :: old_num, new_num

type(curve_plug), pointer :: tempp
type(cv_plug_info) :: plug
type(adss_plug) :: address
integer :: head_mark
logical :: head_switch

ndd(new_num)%status='awake '
ndd(new_num)%n_type=ndd(old_num)%n_type
ndd(new_num)%n_cell=ndd(old_num)%n_cell
ndd(new_num)%cv_rg%total=ndd(old_num)%cv_rg%total
ndd(new_num)%cv_rg%top_idx=ndd(old_num)%cv_rg%top_idx
ndd(new_num)%cv_rg%begin=>ndd(old_num)%cv_rg%begin
ndd(new_num)%cv_rg%eend=>ndd(old_num)%cv_rg%eend
ndd(new_num)%cv_rg%curves_type=ndd(old_num)%cv_rg%curves_type

tempp=>ndd(new_num)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 plug=tempp%plug; address=tempp%address

 address%nd_nmb=new_num; tempp%address=address
 
 select case(plug%end_type)
  case('begin '); cvv(plug%cv_nmb)%begin_end=new_num
  case('end   '); cvv(plug%cv_nmb)%end_end=new_num
  case default; print*, 'Something is wrong!!!'; pause
 end select

 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)
end do

ndd(old_num)%status='asleep'
ndd(old_num)%n_type='   	  '
call clean_up_n_cell(ndd(old_num)%n_cell)
ndd(old_num)%cv_rg%total=0
ndd(old_num)%cv_rg%top_idx=0
nullify(ndd(old_num)%cv_rg%begin)
nullify(ndd(old_num)%cv_rg%eend)
ndd(old_num)%cv_rg%curves_type='      '

end subroutine change_node_number


subroutine evaluation_in_cell(i0,j0,value)
! This subroutine evaluates	the value in a given grid cell.

implicit none
integer, intent(in) :: i0, j0
! The indexes of the given grid cell.
type(state), intent(out) :: value

integer :: nd_nmb
type(node_info) :: n_cell
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: address_c, address_cc, address_ccc
character*6 :: side_c
type(critical_cell), pointer :: temp

value=error_data
select case(ggd_cell(i0,j0)%region)

 case('smth  ')
  value=uu(i0,j0)

 case('node  ')
  do nd_nmb=1,ndn
   if(ndd(nd_nmb)%status.ne.'asleep') then
    n_cell=ndd(nd_nmb)%n_cell
    if(n_cell%x_idx.eq.i0.and.n_cell%y_idx.eq.j0) then
     value=n_cell%or_state
    end if
   end if
  end do
 case('crit  ')
  address_c=ggd_cell(i0,j0)%ccpt(1)%address
  call clean_up_address(address_cc); call clean_up_address(address_ccc)
  side_c=ggd_cell(i0,j0)%ccpt(1)%side
  call visit(address_c,temp)
  g_cell=temp%g_cell; p_cell=temp%p_cell
  value=p_cell%or_state
  select case(side_c)
   case('left  '); address_cc=temp%r_stk; side_c='right '
   case('right '); address_cc=temp%l_stk; side_c='left  '
   case default; call error_message
  end select
  do while(address_cc%cv_nmb.gt.0)
   select case(side_c)
    case('left  '); value=value-p_cell%l_state
    case('right '); value=value-p_cell%r_state
   end select
   call visit(address_cc,temp)
   g_cell=temp%g_cell; p_cell=temp%p_cell
   value=value+p_cell%or_state
   address_ccc=address_cc
   if(temp%l_stk.eq.address_c) then
    side_c='right '; address_cc=temp%r_stk
   else
    if(temp%r_stk.eq.address_c) then
     side_c='left  '; address_cc=temp%l_stk
    else
     call error_message
    end if
   end if		 	 	 
   address_c=address_ccc
  end do         
 case default; call error_message
end select

end subroutine evaluation_in_cell


subroutine clean_node_cells
! This subroutine clean up all the nodes' 'n_cell's.

implicit none

integer :: i, i0, j0

do i=1, ndn
 if(ndd(i)%status.eq.'awake ') then
  i0=ndd(i)%n_cell%x_idx; j0=ndd(i)%n_cell%y_idx
  uu(i0,j0)=error_data
 end if 
end do

end subroutine clean_node_cells    


end module solution