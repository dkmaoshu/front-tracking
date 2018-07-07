module recover_mid_positions
! This module is for recovering middle discontinuity positions 
! from computed cell-averages. In an 'xx'- or 'yy'-type auxiliary
! critical cell, the middle discontinuity position is the 'x'-
! and 'y'-coordinates of the discontinuity position recovered from
! the left, right and ordinary cell-averages. But in an'xy'-type
! critical cell, the middle discontinuity position is 
! discontinuity positions on the two cell-edges.

use solu_comput
! 'solu_com.f90'

!use auxiliary_discontinuity_curves
!! 'auxi_cv.f90'

use pick_discontinuity_positions
! 'Pick_dp.f90'

implicit none

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell, g_cell1, gp, gn
type(phy_info) :: p_cell, pp, pn
type(comp_info) :: c_cell

public  recovery
private temp, g_cell, g_cell1, p_cell, gp, gn, pp, pn, &
        recover_1, recover_2, average, c_cell, pick_positions, tempa, &
		boundary_end


contains


subroutine recovery
! Carry out the recovery procedure for all the 'awake' curves.

implicit none

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call recover_1(acvv(i))
! Carry out the recovery for 'xx'- or 'yy'-type critical cells.

!   call check_list_ac(1,'down  ')

  call average(acvv(i))
! Average the discontinuity positions in every two coupled
! 'xx'- or 'yy'-critical cells.

!   call check_list_ac(1,'down  ')

  call recover_2(acvv(i))

!   call check_list_ac(1,'down  ')

!   call treat_garlic_heads(acvv(i))

 end if
end do

end subroutine recovery


subroutine recover_1(acv)
! Carry out the first recovery in an individual discontinuity curve. This recovers middle
! discontinuity positions for 'xx'- or 'yy'-type critical cells.

implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend%previous%next%address%idx
 end if
end if

do while(associated(temp).and.head_switch) 

 call clean_up_g_cell(g_cell); call clean_up_p_cell(p_cell)
 call clean_up_g_cell(gp); call clean_up_g_cell(gn)
 call clean_up_p_cell(pp); call clean_up_p_cell(pn)

 address=temp%address
 g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
 if(associated(temp%p_nxt_dr)) then
  gp=temp%p_nxt_dr%g_cell; pp=temp%p_nxt_dr%p_cell
 end if
 if(associated(temp%n_nxt_dr)) then
  gn=temp%n_nxt_dr%g_cell; pn=temp%n_nxt_dr%p_cell
 end if

 select case(g_cell%g_type)
  case('xx '); call middle_position_1
  case('yy '); call middle_position_1
! Compute the middle positions for the 'xx'- or 'yy'-type 
! auxiliary critical cells from the computed left, right
! and ordinary cell-averages.

  case('xy ')
  case default; call error_message
 end select

 temp%g_cell=g_cell 

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do


contains


subroutine middle_position_1
! Compute the middle discontinuity positions for 'xx'- or 
! 'yy'-type auxiliary critical cells from the computed
! left, right and ordinary cell-averages.

implicit none
type(state) :: yy, zz, main_jump, left_difference, right_difference
type(state), dimension(-3:3) :: ul, ur
real(8), dimension(2) :: pt1, pt2
real(8) :: xx
real(8), dimension(2) :: normal
integer :: wave_number, i, j

zz=0.5d0*(p_cell%l_state+p_cell%r_state)
yy=p_cell%or_state-zz
main_jump=p_cell%l_state-p_cell%r_state
i=g_cell%x_idx; j=g_cell%y_idx
select case(g_cell%g_type)
 case('xx ')
  call smth_dpr(i,j,temp,ul,'xx ','left  ','critical',.true.,3,'full    ')
  call smth_dpr(i,j,temp,ur,'xx ','right  ','critical',.true.,3,'full    ')
 case('yy ')
  call smth_dpr(i,j,temp,ul,'yy ','left  ','critical',.true.,3,'full    ')
  call smth_dpr(i,j,temp,ur,'yy ','right  ','critical',.true.,3,'full    ')
 case default; call error_message
end select
if(g_cell%point(1).eq.'left  ') then
 left_difference=ul(0)-ul(-1); right_difference=ur(1)-ur(0)
else
 left_difference=ul(1)-ul(0); right_difference=ur(0)-ur(-1)
end if
call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
normal=normal_of_line(pt1,pt2)
wave_number=p_cell%wv_nb

xx=weighted_division(yy,main_jump,left_difference,right_difference, &
                     normal,p_cell%l_state,p_cell%r_state,wave_number)
if(g_cell%g_type.eq.'xx') then
 g_cell%mdis_posi(2)=0.d0
 if(g_cell%point(1).eq.'left') then
  g_cell%mdis_posi(1)=xx
 else
  g_cell%mdis_posi(1)=-xx
 end if
end if

if(g_cell%g_type.eq.'yy') then
 g_cell%mdis_posi(1)=0.d0
 if(g_cell%point(1).eq.'left') then
  g_cell%mdis_posi(2)=xx
 else
  g_cell%mdis_posi(2)=-xx
 end if
end if

end subroutine middle_position_1


end subroutine recover_1


subroutine recover_2(acv)
! Carry out the second recovery in an individual discontinuity curve. This recovers the 
! middle discontinuity  positions for 'xy'-type critical cells.
implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend%previous%next%address%idx
 end if
end if

do while(associated(temp).and.head_switch) 

 call clean_up_g_cell(g_cell); call clean_up_p_cell(p_cell)
 call clean_up_g_cell(gp); call clean_up_g_cell(gn)
! call clean_up_p_cell(pp); call clean_up_p_cell(pn)

 address=temp%address
 g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
 if(associated(temp%p_nxt_dr)) then
  gp=temp%p_nxt_dr%g_cell; pp=temp%p_nxt_dr%p_cell
 end if
 if(associated(temp%n_nxt_dr)) then
  gn=temp%n_nxt_dr%g_cell; pn=temp%n_nxt_dr%p_cell
 end if

 select case(g_cell%g_type)
  case('xx ','yy ')
  case('xy '); call middle_position_2
! compute the middle positions for the 'xy'-type auxiliary
! critical cells.
  case default; call error_message
 end select

 temp%g_cell=g_cell 

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do


contains


subroutine middle_position_2
! compute the middle discontinuity positions for 'xy'-type
! auxiliary critical cells.

implicit none

type(geo_info) :: g_neigh
!type(state) :: d_sides, d_average, d_weight
type(state) :: yy, main_jump, left_difference, right_difference
type(state), dimension(-3:3) :: ul, ur
real(8), dimension(2) :: xi, eta
integer :: i, single
character*6 :: p_side, n_side
character*11 :: view_direction
! The positive and negative sides.
real(8) :: h_neigh, delta, dd, xxx
logical :: if_crossed
real(8), dimension(2) :: normal, normal_1, normal_2, pt1, pt2, nnormal
integer :: wave_number, ii, jj

! real :: xx

! Initialization.
xi=error_data; eta=error_data
call find_single(g_cell,single)
n_side=g_cell%point(single); p_side=side_shift(n_side)

select case(n_side)
 case('left  ')
  main_jump=p_cell%l_state-p_cell%r_state
  yy=p_cell%or_state-p_cell%r_state
 case('right ')
  main_jump=p_cell%r_state-p_cell%l_state
  yy=p_cell%or_state-p_cell%l_state
 case default; call error_message
end select

! Determine whether the discontinuity curve crosses the single grid point.
if_crossed=.false.
call whether_crossed(single,if_crossed)

!if(gp%g_type.eq.'xy '.or.gn%g_type.eq.'xy ') if_crossed=.true.
! Corners in the discontinuity curve should be eliminated.

do i=1,2
! First determines the neighboring critical cell.
 select case(i)
  case(1); g_neigh=gp
  case(2); g_neigh=gn
 end select

 if(g_neigh%g_type.ne.'xx'.and.g_neigh%g_type.ne.'yy') cycle
! The computation is carried out only for 'xx'- or 'yy'-type neighboring critcal cell. 

! Compute the left and right differences.
 ii=g_cell%x_idx; jj=g_cell%y_idx
 select case(g_neigh%g_type)
  case('xx ')
   call smth_dpr(ii,jj,temp,ul,'xx ','left  ','critical',.true.,3,'full    ')
   call smth_dpr(ii,jj,temp,ur,'xx ','right  ','critical',.true.,3,'full    ')
  case('yy ')
   call smth_dpr(ii,jj,temp,ul,'yy ','left  ','critical',.true.,3,'full    ')
   call smth_dpr(ii,jj,temp,ur,'yy ','right  ','critical',.true.,3,'full    ')
  case default; call error_message
 end select
 if(g_cell%point(1).eq.'left  ') then
  left_difference=ul(0)-ul(-1); right_difference=ur(1)-ur(0)
 else
  left_difference=ul(1)-ul(0); right_difference=ur(0)-ur(-1) 
 end if
 wave_number=p_cell%wv_nb

! Compute the normal.
 normal_1=error_data; normal_2=error_data
 if(gp%g_type.ne.'xy '.and.gp%g_type.ne.'non') then
  call pick_point(gp,1,pt1); call pick_point(gp,2,pt2)
  normal_1=normal_of_line(pt1,pt2)
 end if
 if(gn%g_type.ne.'xy '.and.gn%g_type.ne.'non') then
  call pick_point(gn,1,pt1); call pick_point(gn,2,pt2)
  normal_2=normal_of_line(pt1,pt2)
 end if
 if(normal_1(1).gt.0.9d0*error_data) then
  if(normal_2(1).gt.0.9d0*error_data) then
   normal=0.5d0*(normal_1+normal_2)
  else
   normal=normal_1
  end if
 else
  if(normal_2(1).gt.0.9d0*error_data) then
   normal=normal_2
  else
   call error_message
  end if
 end if

! Determine the view-direction, whether it is viewed from the top or right, or from the 
! bottom or left.
 select case(single)
  case(1); view_direction='bottom_left'
  case(2)
   select case(g_cell%edge(i))
    case(1); view_direction='top_right  '
    case(2); view_direction='bottom_left'
    case default; call error_message
   end select
  case(3); view_direction='top_right  '
  case(4)
   select case(g_cell%edge(i))
    case(3); view_direction='bottom_left'
    case(4); view_direction='top_right  '
   end select
  case default; call error_message
 end select
   
! Determine the hight in the neighboring critical cell.
 select case(g_neigh%g_type)
  case('xx '); h_neigh=g_neigh%mdis_posi(1)
  case('yy '); h_neigh=g_neigh%mdis_posi(2)
 end select   
 select case(view_direction)
  case('bottom_left'); h_neigh=h_neigh+0.5d0
  case('top_right  '); h_neigh=0.5d0-h_neigh
 end select

! Now compute the discontinuity positions.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! dd=weighted_division(d_average,d_sides,d_weight)
 dd=weighted_division(yy,main_jump,left_difference,right_difference, &
                      normal,p_cell%l_state,p_cell%r_state,wave_number)
 delta=dabs(dd)
 delta=delta*delta+delta*h_neigh
 if(delta.lt.0.0d0) delta=0.0d0
 delta=dsqrt(delta)

!  dd=dabs(dd)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 if(dd.lt.0.0d0) if_crossed=.true. 
 
 if(.not.if_crossed) then

   dd=dabs(dd)

  eta(i)=-2.0d0*dd+2.0d0*delta
 else
  eta(i)=-2.0d0*dd-2.0d0*delta
 end if
 xi(i)=0.5d0*eta(i)/(h_neigh-eta(i))

 if(if_crossed) then
  if(eta(i).gt.0.0d0) eta(i)=-1.0d-6
  if(xi(i).gt.0.0d0) xi(i)=-1.0d-6
 else
  if(eta(i).lt.0.0d0) eta(i)=1.0d-6
  if(xi(i).lt.0.0d0) xi(i)=1.0d-6
 end if

 select case(view_direction)
  case('bottom_left'); eta(i)=eta(i)-0.5d0
  case('top_right  '); eta(i)=0.5d0-eta(i)
  case default; call error_message
 end select

 select case(g_cell%edge(i))
  case(1,4); xi(i)=xi(i)-0.5d0
  case(2,3); xi(i)=0.5d0-xi(i)
 end select

end do

!g_cell%mdis_posi=g_cell%dis_posi

g_cell%mdis_posi=error_data

do i=1,2
 if(eta(i).gt.0.9d0*error_data) then
  g_cell%mdis_posi(i)=eta(i)
 else
  select case(i)
   case(1)
    if(xi(2).gt.0.9d0*error_data) then
     g_cell%mdis_posi(1)=xi(2)
     if(associated(temp%p_nxt_dr)) then
     if(gp%g_type.ne.'xy ') call error_message
      g_cell%mdis_posi(1)=0.5d0*(gp%mdis_posi(2)+g_cell%mdis_posi(1))
     end if
    end if
   case(2)
    if(xi(1).gt.0.9d0*error_data) g_cell%mdis_posi(2)=xi(1)
  end select
 end if
end do

left_difference=0.0d0; right_difference=0.0d0; xi=error_data

if(g_cell%mdis_posi(1).lt.0.9d0*error_data.or.g_cell%mdis_posi(2).lt.0.9d0*error_data) then

 call normal_4_isolated_xy(g_cell,normal)
 wave_number=p_cell%wv_nb
 dd=weighted_division(yy,main_jump,left_difference,right_difference,normal,p_cell%l_state, &
                      p_cell%r_state,wave_number)
 
 if(dd.lt.0.0d0) if_crossed=.true.
  
 dd=dabs(dd)
  
 nnormal(1)=dabs(normal(1)); nnormal(2)=dabs(normal(2))
 if(nnormal(1).gt.3.0d0*nnormal(2)) then
  nnormal(1)=3.0d0/dsqrt(10.0d0); nnormal(2)=1.0d0/dsqrt(10.0d0)
 end if
 if(3.0d0*nnormal(1).lt.nnormal(2)) then
  nnormal(1)=1.0d0/dsqrt(10.0d0); nnormal(2)=3.0d0/dsqrt(10.0d0)
 end if

 if(g_cell%mdis_posi(1).lt.0.9d0*error_data.and.g_cell%mdis_posi(2).lt.0.9d0*error_data) then
  xxx=dsqrt(2.0d0)/2.0d0
  select case(single)
   case(1)
    if(g_cell%edge(1).eq.1) then
     nnormal=(/xxx,xxx/)
    else
     nnormal=(/-xxx,-xxx/)
    end if
   case(2)
    if(g_cell%dis_posi(1).eq.1) then
     nnormal=(/xxx,-xxx/)
    else
     nnormal=(/-xxx,xxx/)
    end if
   case(3)
    if(g_cell%edge(1).eq.2) then
     nnormal=(/xxx,xxx/)
    else
     nnormal=(/-xxx,-xxx/)
    end if
   case(4)
    if(g_cell%dis_posi(1).eq.1) then
     nnormal=(/xxx,-xxx/)
    else
     nnormal=(/-xxx,xxx/)
    end if
  end select
 end if
   

 xi(1)=dsqrt(2.0d0*dd*nnormal(2)/nnormal(1))
 xi(2)=xi(1)*nnormal(2)/nnormal(1)

 select case(single)

  case(1) 
   do i=1,2
    if(g_cell%mdis_posi(i).lt.0.9d0*error_data) then
     select case(g_cell%edge(i))
      case(1)
       if(g_cell%dis_posi(i).lt.-0.5d0.or.if_crossed) then
        g_cell%mdis_posi(i)=-xi(1)-0.5d0
       else
        g_cell%mdis_posi(i)=xi(1)-0.5d0
       end if
      case(4)
       if(g_cell%dis_posi(i).lt.-0.5d0.or.if_crossed) then
        g_cell%mdis_posi(i)=-xi(2)-0.5d0
       else
        g_cell%mdis_posi(i)=xi(2)-0.5d0
       end if
      case default; call error_message
     end select
    end if
   end do 

  case(2)
   do i=1,2
    if(g_cell%mdis_posi(i).lt.0.9d0*error_data) then 
     select case(g_cell%edge(i)) 
      case(1)
       if(g_cell%dis_posi(i).lt.0.5d0.or..not.if_crossed) then
        g_cell%mdis_posi(i)=-xi(1)+0.5d0
       else
        g_cell%mdis_posi(i)=xi(1)+0.5d0
       end if
      case(2)
       if(g_cell%dis_posi(i).lt.-0.5d0.or.if_crossed) then
        g_cell%mdis_posi(i)=-xi(2)-0.5d0
       else
        g_cell%mdis_posi(i)=xi(2)-0.5d0
       end if
      case default; call error_message
     end select
    end if
   end do 

  case(3) 
   do i=1,2
    if(g_cell%mdis_posi(i).lt.0.9d0*error_data) then 
     select case(g_cell%edge(i)) 
      case(2)
       if(g_cell%dis_posi(i).lt.0.5d0.or..not.if_crossed) then
        g_cell%mdis_posi(i)=-xi(1)+0.5d0
       else
        g_cell%mdis_posi(i)=xi(1)+0.5d0
       end if
      case(3)
       if(g_cell%dis_posi(i).lt.0.5d0.or..not.if_crossed) then
        g_cell%mdis_posi(i)=-xi(2)+0.5d0
       else
        g_cell%mdis_posi(i)=xi(2)+0.5d0
       end if
      case default; call error_message
     end select
    end if
   end do 

  case(4)
   do i=1,2
    if(g_cell%mdis_posi(i).lt.0.9d0*error_data) then 
     select case(g_cell%edge(i)) 
      case(3)
       if(g_cell%dis_posi(i).lt.-0.5d0.or.if_crossed) then
        g_cell%dis_posi(i)=-xi(1)-0.5d0
       else
        g_cell%dis_posi(i)=xi(1)-0.5d0
       end if
       case(4)
       if(g_cell%dis_posi(i).lt.0.5d0.or..not.if_crossed) then
        g_cell%dis_posi(i)=-xi(2)+0.5d0
       else
        g_cell%dis_posi(i)=xi(2)+0.5d0
       end if
      case default; call error_message
     end select
    end if 
   end do 

  case default; call error_message
 end select    

 do i=1, 2
  if(dabs(g_cell%mdis_posi(i)).gt.1.0d0) then
   g_cell%mdis_posi(i)=dsign(1.0d0,g_cell%mdis_posi(i))
  end if
 end do
  
end if

end subroutine middle_position_2


subroutine whether_crossed(single,if_crossed)
! This subroutine determines whether the discontinuity curve should crosse the single
! grid point.

integer, intent(in) :: single
logical, intent(out) :: if_crossed

integer :: i

if_crossed=.false.
select case(single)
 case(1)
  if(g_cell%dis_posi(1).lt.-0.5d0.or.g_cell%dis_posi(2).lt.-0.5d0) then
   if_crossed=.true.
   if(g_cell%dis_posi(1).ge.-0.5d0) g_cell%dis_posi(1)=-0.5000001
   if(g_cell%dis_posi(2).ge.-0.5d0) g_cell%dis_posi(2)=-0.5000001
  end if
 case(2)
  do i=1,2
   if(g_cell%dis_posi(i).gt.0.5d0.and.g_cell%edge(i).eq.1) then
    if_crossed=.true.; exit
   end if
   if(g_cell%dis_posi(i).lt.-0.5d0.and.g_cell%edge(i).eq.2) then
    if_crossed=.true.; exit
   end if
  end do
  if(if_crossed) then
   do i=1,2
    if(g_cell%dis_posi(i).le.0.5d0.and.g_cell%edge(i).eq.1) g_cell%dis_posi(i)=0.5000001
    if(g_cell%dis_posi(i).ge.0.-5d0.and.g_cell%edge(i).eq.2) g_cell%dis_posi(i)=-0.5000001
   end do
  end if
 case(3)
  if(g_cell%dis_posi(1).gt.0.5d0.or.g_cell%dis_posi(2).gt.0.5d0) then
   if_crossed=.true.
   if(g_cell%dis_posi(1).le.0.5d0) g_cell%dis_posi(1)=0.5000001
   if(g_cell%dis_posi(2).le.0.5d0) g_cell%dis_posi(2)=0.5000001
  end if
 case(4)
  do i=1,2
   if(g_cell%dis_posi(i).gt.0.5d0.and.g_cell%edge(i).eq.4) then
    if_crossed=.true.; exit
   end if
   if(g_cell%dis_posi(i).lt.-0.5d0.and.g_cell%edge(i).eq.3) then
    if_crossed=.true.; exit
   end if
  end do
  if(if_crossed) then
   do i=1,2
    if(g_cell%dis_posi(i).le.0.5d0.and.g_cell%edge(i).eq.4) g_cell%dis_posi(i)=0.5000001
    if(g_cell%dis_posi(i).ge.0.-5d0.and.g_cell%edge(i).eq.3) g_cell%dis_posi(i)=-0.5000001
   end do
  end if
 case default; call error_message
end select

end subroutine whether_crossed


end subroutine recover_2


subroutine average(acv)
! Avarage the middle discontinuity positions in coupled auxiliary
! critical cells.

implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch
character*8 :: partner
real(8), dimension(2) :: mposi
integer :: i

if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend%previous%next%address%idx
 end if
end if

do while(associated(temp).and.head_switch) 

 address=temp%address
 g_cell=temp%g_cell; partner=temp%partner
 if(partner.eq.'next    ') then
  temp=>temp%next; g_cell1=temp%g_cell
  
  do i=1, 2
   mposi(i)=0.5d0*(g_cell%mdis_posi(i)+g_cell1%mdis_posi(i))
  end do
  
  if(g_cell%g_type.eq.'xx') then
   if(g_cell%x_idx.gt.g_cell1%x_idx) then
    g_cell%mdis_posi(1)=mposi(1)-0.5d0
    g_cell1%mdis_posi(1)=mposi(1)+0.5d0
   else
    g_cell%mdis_posi(1)=mposi(1)+0.5d0
    g_cell1%mdis_posi(1)=mposi(1)-0.5d0
   end if
  else
   if(g_cell%y_idx.gt.g_cell1%y_idx) then
    g_cell%mdis_posi(2)=mposi(2)-0.5d0
    g_cell1%mdis_posi(2)=mposi(2)+0.5d0
   else
    g_cell%mdis_posi(2)=mposi(2)+0.5d0
    g_cell1%mdis_posi(2)=mposi(2)-0.5d0
   end if
  end if

  temp%previous%g_cell=g_cell; temp%g_cell=g_cell1
 end if
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine average


end module recover_mid_positions