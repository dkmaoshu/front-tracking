module solu_comput
! This module contains the basic variables and utilities used in
! the computation of numerical solution.

use grid_map
! 'grid_map.f90'

use solu_in_smth
! 'solution_smooth.f90

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none

type crit_info_a
 type(adss_info) :: address
 type(geo_info) :: g_cell
 type(phy_info) :: p_cell
 type(state), dimension(0:2) :: l_state, r_state, or_state
 character*6 :: side
end type crit_info_a

type(crit_info_a) :: crit

interface smth_dpr
 module procedure smth_dpr_a
end interface

public  uu, acvv, ggd_cell, smth_dpr, revers_auxi_curve
private pick_go_ord, pick_go_crt, smth_dpr_a, transf, crit, &
        pick_smth 


contains


subroutine smth_dpr_a(i,j,temp,uxy,xy_dir,side,start,if_expand,hg_ord,range_slide)
! This subroutine prepares the data for operations in smooth 
! regions. At grid cells that belong to the computed region, it 
! maintains the original data of the solution, and at grid cells 
! that are out of the computed region, it supplies the extrapolated	
! data from the data of the solution within the region. The 
! variable 'hg_ord' is the highest order of extrapolation +1 to 
! be used.

implicit none

integer, intent(in) :: i, j
! The indexes of the starting grid cell.
type(auxi_crit_cell), pointer :: temp
! Containing the information of the critical cell if starting with
! it.
type(state), dimension(-3:3), intent(out) :: uxy
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

!real(8), dimension(-3:3) :: uxyv

call transf(temp,side)
call pick_smth(i,j,uxy,xy_dir,start)
if(if_expand) call data_fill(uxy,hg_ord,range_slide)

end subroutine smth_dpr_a


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

type(crit_info_a) :: cr
integer :: i0, j0
character*8 :: stoped
integer :: k, kl

do kl=-3,3
 uxy(kl)%value=error_data
 uxy(kl)%gamma=error_data
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
type(crit_info_a), intent(inout) :: cr
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

type(auxi_crit_cell), pointer :: temp_g
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

if(i0.lt.nxll_boundary.or.i0.gt.nxll+nx_boundary-1) stoped='blind'
if(j0.lt.nyll_boundary.or.j0.gt.nyll+ny_boundary-1) stoped='blind'
if(stoped.ne.'blind') then
 if(ggd_cell(i0,j0)%region.eq.'smth') then
  stoped='ordinary'
 else
  if(ggd_cell(i0,j0)%region.eq.'crit') then
   stoped='critical'
   chosen_corner=cycle_4(direction)
   cr%address=ggd_cell(i0,j0)%ccpt(chosen_corner)%address
   call visit(cr%address,temp_g)
   cr%g_cell=temp_g%g_cell; cr%p_cell=temp_g%p_cell
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
type(crit_info_a), intent(inout) :: cr
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

type(auxi_crit_cell), pointer :: temp_g, temp, temp_n, temp_nx
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
call visit(cr%address,temp)
right_ccpt%address=temp%address; right_ccpt%side=cr%side

call clean(also_right_ccpt); nullify(temp_n)
select case(cr%side)
 case('left  ')
  if(temp%l_stk%cv_nmb.gt.0) call visit(temp%l_stk,temp_n)
 case('right ')
  if(temp%r_stk%cv_nmb.gt.0) call visit(temp%r_stk,temp_n)
end select
if(associated(temp_n)) then
 also_right_ccpt%address=temp_n%address
 if(temp_n%l_stk.eq.temp%address) then
  also_right_ccpt%side='left  '
 else
  if(temp_n%r_stk.eq.temp%address) then
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

if(i0.lt.nxll_boundary.or.i0.gt.nxll+nx_boundary-1) stoped='blind'
if(j0.lt.nyll_boundary.or.j0.gt.nyll+ny_boundary-1) stoped='blind'
if(stoped.ne.'blind') then
 if(ggd_cell(i0,j0)%region.eq.'smth') then
  stoped='ordinary'
 else
  if(ggd_cell(i0,j0)%region.eq.'crit') then
   stoped='critical'
   cr%address=ggd_cell(i0,j0)%ccpt(chosen_corner)%address
   call visit(cr%address,temp_g)
   cr%g_cell=temp_g%g_cell; cr%p_cell=temp_g%p_cell
!   cr%side=cr%g_cell%point(chosen_corner)
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
 
 select case(direction)
  case(1); i0=i0+1; dir='xx '
  case(2); j0=j0+1; dir='yy '
  case(3); i0=i0-1; dir='xx '
  case(4); j0=j0-1; dir='yy '
 end select
  
 chosen_corner=error_index
 call visit(cr%address,temp); nullify(temp_nx)
 if(temp%partner.ne.'single  ') then
  
  select case(temp%partner)
   case('previous'); temp_nx=>temp%previous
   case('next    '); temp_nx=>temp%next
   case default; call error_message
  end select
  if(associated(temp_nx)) then
   gg_cl=temp_nx%g_cell; pp_cl=temp_nx%p_cell
   if(i0.eq.gg_cl%x_idx.and.j0.eq.gg_cl%y_idx) then!
	cr%g_cell=gg_cl; cr%p_cell=pp_cl
	stoped='critical'
   end if
  endif
 end if   

 if(stoped.eq.'blind   ') then
   
  do k=1,2
   if(cr%g_cell%point(corners(k)).eq.cr%side) then
    chosen_corner=corners(k); exit
   end if
  end do
   
  if(chosen_corner.gt.error_index) then
   
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
 end if  
endif   

end subroutine pick_go_crt


subroutine transf(temp,side)

implicit none

type(auxi_crit_cell), pointer :: temp
character*6, intent(in) :: side

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
 crit%g_cell=temp%g_cell; crit%p_cell=temp%p_cell
 crit%side=side
! crit%l_state%value=error_data; crit%r_state%value=error_data
! crit%or_state%value=error_data
end if

end subroutine transf


subroutine revers_auxi_curve(curve_number,if_revers_g,if_revers_p, &
                             if_revers_st)
! This subroutine reverses an auxiliary discontinuity curve list.

implicit none
integer, intent(in) :: curve_number
! The number of the curve to be reversed.
character*3, intent(in) :: if_revers_g
! Indicate whether the geometrical information in each critical
! cell should be reversed.
character*3, intent(in) :: if_revers_p
! Indicate whether the physical information in each critical cell
! should be reversed.
character*3, intent(in) :: if_revers_st
! Indicate whether the left- and right-stack relation should be
! reversed.

type(auxi_crit_cell), pointer :: temp, temp_tp, current
type(geo_info) :: g_cell, g_temp
type(phy_info) :: p_cell, p_temp
type(adss_info) :: address
integer :: head_mark, i, i0, j0, acurve_temp
logical :: head_switch

!  call check_list_ac(1,'down  ')

do i=1,cvn
 if(acvv(i)%status.ne.'awake '.and.acvv(i)%status.ne.'yawn  ')  exit
end do
acurve_temp=i

acvv(acurve_temp)%status='awake '
acvv(acurve_temp)%cv_type=acvv(curve_number)%cv_type
acvv(acurve_temp)%begin_end=acvv(curve_number)%end_end
acvv(acurve_temp)%end_end=acvv(curve_number)%begin_end
acvv(acurve_temp)%wave=wave_revers(acvv(curve_number)%wave)
acvv(acurve_temp)%total=0
call creat_acv(acurve_temp)

current=>acvv(acurve_temp)%eend
if(associated(acvv(curve_number)%begin%next)) then
 head_mark=acvv(curve_number)%begin%next%address%idx
 head_switch=.true.
 temp=>acvv(curve_number)%begin%next
 do while(head_switch)
  g_cell=temp%g_cell; p_cell=temp%p_cell; address=temp%address
  call clean_up_g_cell(g_temp); call clean_up_p_cell(p_temp)
! Firstly, revers the geometrical information if required.
  if(if_revers_g.eq.'yes') then
   g_temp%g_type=g_cell%g_type
   g_temp%x_idx=g_cell%x_idx; g_temp%y_idx=g_cell%y_idx
   g_temp%dis_posi(1)=g_cell%dis_posi(2)
   g_temp%dis_posi(2)=g_cell%dis_posi(1)
   g_temp%edge(1)=g_cell%edge(2); g_temp%edge(2)=g_cell%edge(1)
   do i=1,4
    select case(g_cell%point(i))
     case('left  '); g_temp%point(i)='right '
     case('right '); g_temp%point(i)='left  '
     case default; call error_message
    end select
   end do   	 
  else
   g_temp=g_cell
  end if
! Secondly, revers the physical information.
  if(if_revers_p.eq.'yes') then
   p_temp%l_state=p_cell%r_state
   p_temp%r_state=p_cell%l_state
   p_temp%or_state=p_cell%or_state
   p_temp%wv_nb=acvv(acurve_temp)%wave
  else
   p_temp=p_cell
  end if
! Thirdly, update the grid map.
  i0=g_temp%x_idx; j0=g_temp%y_idx
  do i=1,4
   if(ggd_cell(i0,j0)%ccpt(i)%address.eq.address) then
    ggd_cell(i0,j0)%ccpt(i)%side=g_temp%point(i)
   end if
  end do
! Finally, put the critical cell on the the new cruve list and 
! update the neighboring relation.
  allocate(temp_tp)
  temp_tp%address%cv_nmb=curve_number
  temp_tp%g_cell=g_temp; temp_tp%p_cell=p_temp
  select case(if_revers_st) 
   case('yes')
    temp_tp%l_stk=temp%r_stk; temp_tp%r_stk=temp%l_stk 
   case('no ')
    temp_tp%l_stk=temp%l_stk; temp_tp%r_stk=temp%r_stk
   case default; call error_message
  end select
  call insert(current,temp_tp,'before')
  
  temp=>temp%next
  if(associated(temp)) then
   head_switch=(temp%address%idx.ne.head_mark)
  else
   exit
  end if
 end do  
end if

call deletee(acvv(curve_number))
acvv(curve_number)=acvv(acurve_temp)

acvv(acurve_temp)%status='asleep'
acvv(acurve_temp)%cv_type='   '
acvv(acurve_temp)%begin_end=0
acvv(acurve_temp)%end_end=0
acvv(acurve_temp)%total=0
acvv(acurve_temp)%wave=0
nullify(acvv(acurve_temp)%begin); nullify(acvv(acurve_temp)%eend)

end subroutine revers_auxi_curve


subroutine find_original(temp,g_orig)

implicit none
type(auxi_crit_cell), pointer :: temp
type(geo_info), intent(out) :: g_orig

type(geo_info) :: g_cell, g_partner
integer :: num

g_cell=temp%g_cell; g_orig=g_cell; call clean_up_g_cell(g_partner)
select case(temp%partner)
 case('single  '); g_orig=g_cell
 case('previous','next    ')
  select case(temp%partner)
   case('previous')
    num=1; g_partner=temp%previous%g_cell
   case('next    ')
    num=2; g_partner=temp%next%g_cell
  end select
  select case(g_cell%g_type)
   case('xx ')
    if(g_partner%x_idx.gt.g_cell%x_idx) then
     select case(g_cell%edge(num))
      case(1); g_orig%point(2)=g_orig%point(1)
      case(3); g_orig%point(3)=g_orig%point(1) 
      case default; call error_message
     end select
    else
     select case(g_cell%edge(num))
      case(1); g_orig%point(1)=g_orig%point(3)
      case(3); g_orig%point(4)=g_orig%point(3) 
      case default; call error_message
     end select
    end if
   case('yy ')
    if(g_partner%y_idx.gt.g_cell%y_idx) then
     select case(g_cell%edge(num))
      case(2); g_orig%point(3)=g_orig%point(1)
      case(4); g_orig%point(4)=g_orig%point(1) 
      case default; call error_message
     end select
    else
     select case(g_cell%edge(num))
      case(2); g_orig%point(2)=g_orig%point(3)
      case(4); g_orig%point(1)=g_orig%point(3) 
      case default; call error_message
     end select
    end if
  end select
 case default; call error_message
end select

end subroutine find_original


end module solu_comput
