module clean_up_2_4_contact
! This module is for the clean_up procedure of stacked discontinuity cells 
! of contact discontinuities.

use discontinuity_curves
! 'discontinuity_curve.f90'

implicit none


contains


subroutine clean_up_2
! Implement the clean_up.

implicit none

integer :: i

! call check_list_c(7,'down  ')

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call clean_memo(i)
 end if
end do  

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call cleaning(i)
 end if
end do

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call return_data(i)
 end if
end do


contains


subroutine clean_memo(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp
type(adss_info) :: address
integer :: head_mark
logical :: head_switch

nullify(temp)
if(associated(cvv(cv_nmb)%begin%next)) then
 temp=>cvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 address=temp%address
 
 temp%operation_memo=.false.
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine clean_memo


subroutine cleaning(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp, temp_n
type(adss_info) :: address
integer :: head_mark
logical :: head_switch

type(geo_info) :: g_cell, g_cell_n
type(phy_info) :: p_cell, p_cell_n
type(state) :: l_store, r_store, l_store_n, r_store_n, ordin_state
type(state) :: left_phy_state, right_phy_state, left_phy_state_n, &
               right_phy_state_n, ordin_phy_state
real(8) :: l_x_velo, r_x_velo, l_y_velo, r_y_velo, l_x_velo_n, r_x_velo_n, &
           l_y_velo_n, r_y_velo_n, or_x_velo, or_y_velo
real(8) :: l_pressure, r_pressure, l_pressure_n, r_pressure_n, &
           or_pressure, ttt
character*6 :: side, side_n

nullify(temp)
if(associated(cvv(cv_nmb)%begin%next)) then
 temp=>cvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell
 address=temp%address

 if(temp%l_stk.ne.adss_info(0,0).or.temp%r_stk.ne.adss_info(0,0)) then
! Update the left and right states only for stacked discontinuity cells.
  
  nullify(temp_n)
  if(temp%l_stk.ne.adss_info(0,0)) then
   call visit(temp%l_stk,temp_n); side='left  ' 
  else
   call visit(temp%r_stk,temp_n); side='right '
  end if

  if(temp%l_stk.ne.adss_info(0,0).and.temp%r_stk.ne.adss_info(0,0)) then
   call error_message
  end if
  if(temp_n%l_stk.ne.adss_info(0,0).and.temp_n%r_stk.ne.adss_info(0,0)) then
   call error_message
  end if
! Two-side stacked case has not been treated in the present code, and it
! will be done the future version.
   
  if(temp%p_cell%wv_nb.ne.2.or.temp_n%p_cell%wv_nb.ne.2) goto 10
! Only stacked contact discontinuity cells are treated.
   
  if(temp_n%l_stk.eq.temp%address) then
   side_n='left  '
  else
   if(temp_n%r_stk.eq.temp%address) then
    side_n='right '
   else
    call error_message
   end if
  end if

  g_cell_n=temp_n%g_cell; p_cell_n=temp_n%p_cell

  ordin_state=p_cell%or_state+p_cell_n%or_state
  if(side.eq.'left  ') then
   ordin_state=ordin_state-p_cell%l_state
  else
   ordin_state=ordin_state-p_cell%r_state
  end if
  
  call find_physical_state(p_cell%l_state,(/1.0d0,0.0d0/),left_phy_state)
  l_x_velo=left_phy_state%value(2)
  l_y_velo=left_phy_state%value(3)
  l_pressure=left_phy_state%value(4)
  call find_physical_state(p_cell%r_state,(/1.0d0,0.0d0/),right_phy_state)
  r_x_velo=right_phy_state%value(2)
  r_y_velo=right_phy_state%value(3)
  r_pressure=right_phy_state%value(4)

  call find_physical_state(p_cell_n%l_state,(/1.0d0,0.0d0/),left_phy_state_n)
  l_x_velo_n=left_phy_state_n%value(2)
  l_y_velo_n=left_phy_state_n%value(3)
  l_pressure_n=left_phy_state_n%value(4)
  call find_physical_state(p_cell_n%r_state,(/1.0d0,0.0d0/),right_phy_state_n)
  r_x_velo_n=right_phy_state_n%value(2)
  r_y_velo_n=right_phy_state_n%value(3)
  r_pressure_n=right_phy_state_n%value(4)

  call find_physical_state(ordin_state,(/1.0d0,0.0d0/),ordin_phy_state)
  or_x_velo=ordin_phy_state%value(2)
  or_y_velo=ordin_phy_state%value(3)
  or_pressure=ordin_phy_state%value(4) 
   
  ttt=or_x_velo-l_x_velo
  l_x_velo=l_x_velo+ttt
  ttt=or_x_velo-r_x_velo
  r_x_velo=r_x_velo+ttt
  ttt=or_x_velo-l_x_velo_n
  l_x_velo_n=l_x_velo+ttt
  ttt=or_x_velo-r_x_velo_n
  r_x_velo_n=r_x_velo+ttt

  ttt=or_y_velo-l_y_velo
  l_y_velo=l_y_velo+ttt
  ttt=or_y_velo-r_y_velo
  r_y_velo=r_y_velo+ttt
  ttt=or_y_velo-l_y_velo_n
  l_y_velo_n=l_y_velo+ttt
  ttt=or_y_velo-r_y_velo_n
  r_y_velo_n=r_y_velo+ttt
    
  ttt=or_pressure-l_pressure
  l_pressure=l_pressure+ttt
  ttt=or_pressure-r_pressure
  r_pressure=r_pressure+ttt
  ttt=or_pressure-l_pressure_n
  l_pressure_n=l_pressure+ttt
  ttt=or_pressure-r_pressure_n
  r_pressure_n=r_pressure+ttt
  
  left_phy_state%value(2)=l_x_velo
  left_phy_state%value(3)=l_y_velo
  left_phy_state%value(4)=l_pressure
  right_phy_state%value(2)=r_x_velo
  right_phy_state%value(3)=r_y_velo
  right_phy_state%value(4)=r_pressure
  left_phy_state_n%value(2)=l_x_velo_n
  left_phy_state_n%value(3)=l_y_velo_n
  left_phy_state_n%value(4)=l_pressure_n
  right_phy_state_n%value(2)=r_x_velo_n
  right_phy_state_n%value(3)=r_y_velo_n
  right_phy_state_n%value(4)=r_pressure_n
  
  call find_conservative_state(left_phy_state,(/1.0d0,0.0d0/),l_store)
  call find_conservative_state(right_phy_state,(/1.0d0,0.0d0/),r_store)
  call find_conservative_state(left_phy_state_n,(/1.0d0,0.0d0/),l_store_n)
  call find_conservative_state(right_phy_state_n,(/1.0d0,0.0d0/),r_store_n)
    
 end if
  
 temp%l_store=l_store; temp%r_store=r_store
 temp_n%l_store=l_store_n; temp_n%r_store=r_store_n

 10  temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine cleaning


end subroutine clean_up_2


end module clean_up_2_4_contact