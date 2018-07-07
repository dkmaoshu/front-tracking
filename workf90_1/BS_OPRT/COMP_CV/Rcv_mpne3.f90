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
private temp, g_cell, g_cell1, gp, gn, p_cell, recover_1, recover_2, average, c_cell


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
! call clean_up_g_cell(gp); call clean_up_g_cell(gn)
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

 call clean_up_g_cell(gn); call clean_up_p_cell(pn)

 address=temp%address
 g_cell=temp%g_cell; p_cell=temp%p_cell 
 gn=temp%next%g_cell; pn=temp%next%p_cell 
 if(g_cell%g_type.eq.'xy ') then
  call middle_position_20 
! Treatment of neighoring 'xy'-type discontinuity cells when one of them has 'cross' situation.
  temp%g_cell=g_cell; temp%p_cell=p_cell
  temp%next%g_cell=gn; temp%next%p_cell=pn
 end if

! compute the middle positions for the 'xy'-type auxiliary critical cells.

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

! call check_list_ac(1,'down  ')

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

! call clean_up_g_cell(g_cell); call clean_up_p_cell(p_cell)
 call clean_up_g_cell(gp);  call clean_up_g_cell(gn); call clean_up_p_cell(pn)

 address=temp%address
 g_cell=temp%g_cell; p_cell=temp%p_cell !; c_cell=temp%c_cell
 if(associated(temp%p_nxt_dr)) then
  gp=temp%p_nxt_dr%g_cell !; pp=temp%p_nxt_dr%p_cell
 end if
 if(associated(temp%n_nxt_dr)) then
  gn=temp%n_nxt_dr%g_cell !; pn=temp%n_nxt_dr%p_cell
 end if

 if(g_cell%g_type.eq.'xy ') then
  call middle_position_2
 end if
! compute the middle positions for the 'xy'-type auxiliary
! critical cells.

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

! call check_list_ac(1,'down  ')

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
 g_cell=temp%g_cell 

 if(g_cell%g_type.eq.'xy ') then
  call middle_position_21; temp%g_cell=g_cell
 end if

! compute the middle positions for the 'xy'-type auxiliary critical cells.

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do


contains


subroutine middle_position_2
! Compute a point on the discontinuity segment for an 'xy'-type auxiliary critical cell. The point is stored on 'temp%memo_point'.

implicit none

!type(geo_info) :: g_neigh
!type(state) :: d_sides, d_average, d_weight
type(state) :: yy, main_jump, left_difference, right_difference
type(state), dimension(-3:3) :: ul, ur
!real(8), dimension(2) :: xi, eta
integer :: i, single
character*6 :: p_side, n_side
!character*11 :: view_direction
! The positive and negative sides.
real(8) :: dd
!logical :: if_crossed
real(8), dimension(2) :: normal, normal_0, normal_1, normal_2, pt1, pt2 
real(8) :: length_0, length_1, length_2, length
integer :: wave_number, ii, jj

!  real :: xx

! Initialization.
!xi=error_data; eta=error_data
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

!! Determine whether the discontinuity curve crosses the single grid point.
!if_crossed=.false.
!call whether_crossed(single,if_crossed)
wave_number=p_cell%wv_nb
 
! Compute the normal.
normal_0=0.0d0; normal_1=0.0d0; normal_2=0.0d0
! if(gp%g_type.ne.'xy '.and.gp%g_type.ne.'non') then
call pick_point(gp,1,pt1); call pick_point(gp,2,pt2); length_1=length_of_segment(pt1,pt2)
if(length_1.gt.1.0d-6) normal_1=normal_of_line(pt1,pt2)
call pick_point(gn,1,pt1); call pick_point(gn,2,pt2); length_2=length_of_segment(pt1,pt2)
if(length_2.gt.1.0d-6) normal_2=normal_of_line(pt1,pt2)
call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2); length_0=length_of_segment(pt1,pt2)
if(length_0.gt.1.0d-6) normal_0=normal_of_line(pt1,pt2)
if(length_0.lt.1.0d-6.and.length_1.lt.1.0d-6.and.length_2.lt.1.0d-6) call error_message
normal=(length_0*normal_0+length_1*normal_1+length_2*normal_2)/(length_0+length_1+length_2)
length=dsqrt(normal(1)*normal(1)+normal(2)*normal(2))
normal=normal/length
  
! Now compute the middle point from the conservation ratio. This is the real one that is on the discontinuity segment.
dd=weighted_division(yy,main_jump,left_difference,right_difference,normal,p_cell%l_state,p_cell%r_state,wave_number)
dd=dsign(1.0d0,dd)*dsqrt(2.0d0*dabs(dd*normal(1)*normal(2)))

select case(single)
 case(1)
  if(g_cell%point(single).eq.'left  ') then
   if(normal(1).lt.0.0d0.and.normal(2).lt.0.0d0) call error_message 
   if(normal(1).lt.0.0d0) then
    normal(1)=0.0d0; normal(2)=1.0d0; print*, '1_111'
   end if
   if(normal(2).lt.0.0d0) then
    normal(2)=0.0d0; normal(1)=1.0d0; print*, '1_222'
   end if	 	  
   temp%memo_point=(/-0.5d0,-0.5d0/)+dd*normal
  else
   if(g_cell%point(single).ne.'right ') call error_message
   if(normal(1).gt.0.0d0.and.normal(2).gt.0.0d0) call error_message 
   if(normal(1).gt.0.0d0) then
    normal(1)=0.0d0; normal(2)=-1.0d0; print*, '1_333'
   end if
   if(normal(2).gt.0.0d0) then
    normal(2)=0.0d0; normal(1)=-1.0d0; print*, '1_444'
   end if
   temp%memo_point=(/-0.5d0,-0.5d0/)-dd*normal
  end if
 case(2) 
  if(g_cell%point(single).eq.'left  ') then
   if(normal(1).gt.0.0d0.and.normal(2).lt.0.0d0) call error_message
   if(normal(1).gt.0.0d0) then
    normal(1)=0.0d0; normal(2)=1.0d0; print*, '2_111'
   end if
   if(normal(2).lt.0.0d0) then
    normal(2)=0.0d0; normal(1)=-1.0d0; print*, '2_222'
   end if
   temp%memo_point=(/0.5d0,-0.5d0/)+dd*normal
  else 
   if(g_cell%point(single).ne.'right ') call error_message
   if(normal(1).lt.0.0d0.and.normal(2).gt.0.0d0) call error_message
   if(normal(1).lt.0.0d0) then
    normal(1)=0.0d0; normal(2)=-1.0d0; print*, '2_333'
   end if
   if(normal(2).gt.0.0d0) then
    normal(2)=0.0d0; normal(1)=1.0d0; print*, '2_444'
   end if
   temp%memo_point=(/0.5d0,-0.5d0/)-dd*normal
  end if
 case(3) 
  if(g_cell%point(single).eq.'left  ') then
   if(normal(1).gt.0.0d0.and.normal(2).gt.0.0d0) call error_message
   if(normal(1).gt.0.0d0) then
    normal(1)=0.0d0; normal(2)=-1.0d0; print*, '3_111'
   end if
   if(normal(2).gt.0.0d0) then
    normal(2)=0.0d0; normal(1)=-1.0d0; print*, '3_222'
   end if
   temp%memo_point=(/0.5d0,0.5d0/)+dd*normal
  else
   if(g_cell%point(single).ne.'right ') call error_message
   if(normal(1).lt.0.0d0.and.normal(2).lt.0.0d0) call error_message
   if(normal(1).lt.0.0d0) then
    normal(1)=0.0d0; normal(2)=1.0d0; print*, '3_333'
   end if
   if(normal(2).lt.0.0d0) then
    normal(2)=0.0d0; normal(1)=1.0d0; print*, '3_444'
   end if
   temp%memo_point=(/0.5d0,0.5d0/)-dd*normal
  end if
 case(4)
  if(g_cell%point(single).eq.'left  ') then
   if(normal(1).lt.0.0d0.and.normal(2).gt.0.0d0) call error_message
   if(normal(1).lt.0.0d0) then
    normal(1)=0.0d0; normal(2)=-1.0d0; print*, '4_111'
   end if
   if(normal(2).gt.0.0d0) then
    normal(2)=0.0d0; normal(1)=1.0d0; print*, '4_222'
   end if
   temp%memo_point=(/-0.5d0,0.5d0/)+dd*normal
  else 
   if(g_cell%point(single).ne.'right ') call error_message
   if(normal(1).gt.0.0d0.and.normal(2).lt.0.0d0) call error_message
   if(normal(1).gt.0.0d0) then
    normal(1)=0.0d0; normal(2)=1.0d0; print*, '4_333'
   end if
   if(normal(2).lt.0.0d0) then
    normal(2)=0.0d0; normal(1)=-1.0d0; print*, '4_444'
   end if	      	    
   temp%memo_point=(/-0.5d0,0.5d0/)-dd*normal
  end if
 case default; call error_message 
end select

end subroutine middle_position_2


subroutine middle_position_21
! compute the middle discontinuity positions for 'xy'-type auxiliary critical cells.

implicit none

type(geo_info) :: g_neigh
integer :: l, single
real(8), dimension(2,2) :: pts1, pts2
real(8), dimension(2) :: pt
character*10 :: err_mesg 

do l=1, 2

! Choose the discontinuity edge segment.
 select case(g_cell%edge(l))
  case(1); pts1(1,1)=-0.5d0; pts1(1,2)=-0.5d0; pts1(2,1)=0.5d0; pts1(2,2)=-0.5d0
  case(2); pts1(1,1)=0.5d0; pts1(1,2)=-0.5d0; pts1(2,1)=0.5d0; pts1(2,2)=0.5d0
  case(3); pts1(1,1)=0.5d0; pts1(1,2)=0.5d0; pts1(2,1)=-0.5d0; pts1(2,2)=0.5d0
  case(4); pts1(1,1)=-0.5d0; pts1(1,2)=0.5d0; pts1(2,1)=-0.5d0; pts1(2,2)=-0.5d0
 end select

! Choose the middle points in the current 'xy'-cell and the neighboring cell.
 pts2(1,:)=temp%memo_point
  
 if(l.eq.1) then
  g_neigh=temp%p_nxt_dr%g_cell
 else 
  g_neigh=temp%n_nxt_dr%g_cell
 end if
 select case(g_neigh%g_type)
  case('xx ','yy '); pts2(2,:)=g_neigh%mdis_posi
  case('xy ') 
   if(l.eq.1) pts2(2,:)=temp%p_nxt_dr%memo_point
   if(l.eq.2) pts2(2,:)=temp%n_nxt_dr%memo_point
  case default; call error_message
 end select
 pts2(2,:)=pts2(2,:)-(/dfloat(g_cell%x_idx-g_neigh%x_idx),dfloat(g_cell%y_idx-g_neigh%y_idx)/)
 call ints_of_lines(pts1,pts2,pt,err_mesg)
 if(err_mesg.eq.'short     ') pt=0.5d0*(pts2(1,:)+pts2(2,:))
 if(err_mesg.eq.'parallel  ') call error_message
  
 select case(g_cell%edge(l))
  case(1,3); g_cell%mdis_posi(l)=pt(1)
  case(2,4); g_cell%mdis_posi(l)=pt(2)
  case default; call error_message
 end select
  
end do

end subroutine middle_position_21


subroutine middle_position_20

implicit none

integer :: corner1, corner2
character*3 :: direction
character*6 :: side
real(8) :: dd1, dd2, ddd1, ddd2
type(state) :: state1, state2

if(gn%g_type.ne.'xy ') return
if(temp%next%next%g_cell%g_type.eq.'xy ') call error_message

call find_single(g_cell,corner1); call find_single(gn,corner2)
if(g_cell%x_idx.eq.gn%x_idx) then
 direction='yy '
else
 if(g_cell%y_idx.eq.gn%y_idx) then
  direction='xx '
 else
  call error_message
 end if
end if 
if(point_shift(corner1,direction).ne.corner2) call error_message  
if(g_cell%point(corner1).ne.gn%point(corner2)) call error_message
! Above is the geometrical check: the corners of the two 'xy'-type discontinuity cells are on the same grid point and with the same side.

side=g_cell%point(corner1)
select case(side)
 case('left  ')
  dd1=p_cell%or_state%value(1)-p_cell%r_state%value(1); dd2=pn%or_state%value(1)-pn%r_state%value(1)
  state1=p_cell%r_state; state2=pn%r_state
 case('right ')
  dd1=p_cell%or_state%value(1)-p_cell%l_state%value(1); dd2=pn%or_state%value(1)-pn%l_state%value(1)
  state1=p_cell%l_state; state2=pn%l_state
 case default; call error_message
end select
if(dd1*dd2.gt.0.0d0) return ! If either discontinuity cell is not 'cross' or both are 'cross', there is no treatmnet.

print*, 'Neighboring xy-type', g_cell%x_idx, g_cell%y_idx, gn%x_idx, gn%y_idx

if((p_cell%l_state%value(1)-p_cell%r_state%value(1))*(pn%l_state%value(1)-pn%r_state%value(1)).le.0.0d0) call error_message
if(p_cell%l_state%value(1)-p_cell%r_state%value(1).gt.0.0d0) then
 if(dd1.gt.0.0d0) then
  p_cell%or_state%value(1)=p_cell%or_state%value(1)-dd1-1.0d-6
  temp%previous%p_cell%or_state%value(1)=temp%previous%p_cell%or_state%value(1)+dd1+1.0d-6
  if(temp%previous%partner.eq.'previous') temp%previous%previous%p_cell%or_state%value(1)=temp%previous%previous%p_cell%or_state%value(1)+dd1+1.0d-6
 else
  pn%or_state%value(1)=pn%or_state%value(1)-dd2-1.0d-6
  temp%next%p_cell%or_state%value(1)=temp%next%p_cell%or_state%value(1)+dd2+1.0d-6
  if(temp%next%partner.eq.'next    ') temp%next%next%p_cell%or_state%value(1)=temp%next%next%p_cell%or_state%value(1)+dd2+1.0d-6
 end if
else
 if(dd1.lt.0.0d0) then
  p_cell%or_state%value(1)=p_cell%or_state%value(1)-dd1+1.0d-6
  temp%previous%p_cell%or_state%value(1)=temp%previous%p_cell%or_state%value(1)+dd1-1.0d-6
  if(temp%previous%partner.eq.'previous') temp%previous%previous%p_cell%or_state%value(1)=temp%previous%previous%p_cell%or_state%value(1)+dd1-1.0d-6
 else
  pn%or_state%value(1)=pn%or_state%value(1)-dd2+1.0d-6
  temp%next%p_cell%or_state%value(1)=temp%next%p_cell%or_state%value(1)+dd2-1.0d-6
  if(temp%next%partner.eq.'next    ') temp%next%next%p_cell%or_state%value(1)=temp%next%next%p_cell%or_state%value(1)+dd2-1.0d-6
 end if
end if
! Transfer the mass. 

end subroutine middle_position_20


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
real(8) :: l_area, r_area, length, l_area1, r_area1, length1
real(8) :: alpha, beta

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
  
  call side_area_length_partnered(g_cell,g_cell1,l_area,r_area,length)
  call side_area_length_partnered(g_cell1,g_cell,l_area1,r_area1,length1)

! The two mid_positions are weighted-averaged with the weights as the lengthes of 
! discontinuity curve segments in the two partnered discontinuity cells.  

  alpha=length/(length+length1); beta=1.0d0-alpha
   
  do i=1, 2
   mposi(i)=alpha*g_cell%mdis_posi(i)+beta*g_cell1%mdis_posi(i)
  end do
  
  if(g_cell%g_type.eq.'xx') then
   if(g_cell%x_idx.gt.g_cell1%x_idx) then
    g_cell%mdis_posi(1)=mposi(1)-beta
    g_cell1%mdis_posi(1)=mposi(1)+alpha
   else
    g_cell%mdis_posi(1)=mposi(1)+beta
    g_cell1%mdis_posi(1)=mposi(1)-alpha
   end if
  else
   if(g_cell%y_idx.gt.g_cell1%y_idx) then
    g_cell%mdis_posi(2)=mposi(2)-beta
    g_cell1%mdis_posi(2)=mposi(2)+alpha
   else
    g_cell%mdis_posi(2)=mposi(2)+beta
    g_cell1%mdis_posi(2)=mposi(2)-alpha
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