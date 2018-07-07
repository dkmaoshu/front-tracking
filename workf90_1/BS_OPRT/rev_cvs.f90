module revers_discontinuity_curves

use solution
! 'solution.f90

public  revers_curves
private revers_curve


contains


subroutine revers_curves

implicit none

integer i

allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ndd(ndn))

call give_cv
call give_gcl
call give_nd

! call check_list_c(1,'down  ')

do i=1,cvn
 if(cvv(i)%status.eq.'awake ') then
  call revers_curve(i)
 end if
end do

! call check_list_c(1,'down  ')

call get_cv
call get_gcl
call get_nd

deallocate(cvv); deallocate(ggd_cell); deallocate(ndd)

end subroutine revers_curves


subroutine revers_curve(curve_number)
! This subroutine reverses a discontinuity curve list.

implicit none
integer, intent(in) :: curve_number
! The number of the curve to be reversed.

type(critical_cell), pointer :: temp, temp_tp, current
type(geo_info) :: g_cell, g_temp
type(phy_info) :: p_cell, p_temp
type(adss_info) :: address
integer :: head_mark, i, i0, j0, curve_temp
logical :: head_switch

!  call check_list_ac(1,'down  ')

do i=1,cvn
 if(cvv(i)%status.ne.'awake ') exit
end do
curve_temp=i

cvv(curve_temp)%status='awake '
cvv(curve_temp)%cv_type=cvv(curve_number)%cv_type
cvv(curve_temp)%begin_end=cvv(curve_number)%end_end
cvv(curve_temp)%end_end=cvv(curve_number)%begin_end
cvv(curve_temp)%total=cvv(curve_number)%total
cvv(curve_temp)%wave=wave_revers(cvv(curve_number)%wave)
call creat_cv(curve_temp)

current=>cvv(curve_temp)%eend
if(associated(cvv(curve_number)%begin%next)) then
 head_mark=cvv(curve_number)%begin%next%address%idx
 head_switch=.true.
 temp=>cvv(curve_number)%begin%next
 do while(head_switch)
  g_cell=temp%g_cell; p_cell=temp%p_cell; address=temp%address
  call clean_up_g_cell(g_temp); call clean_up_p_cell(p_temp)
! Firstly, revers the geometrical information if required.
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
! Secondly, revers the physical information.
  p_temp%l_state=p_cell%r_state
  p_temp%r_state=p_cell%l_state
  p_temp%or_state=p_cell%or_state
  p_temp%wv_nb=cvv(curve_temp)%wave
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
  temp_tp%address=address
  temp_tp%g_cell=g_temp; temp_tp%p_cell=p_temp
  temp_tp%l_stk=temp%r_stk; temp_tp%r_stk=temp%l_stk 

!   call check_list_c(1,'down  ')
      
  call insertt
  current=>temp_tp
   
!   call check_list_c(1,'down  ')
   
  temp=>temp%next
  if(associated(temp)) then
   head_switch=(temp%address%idx.ne.head_mark)
  else
   exit
  end if
 end do  
end if

if(cvv(curve_temp)%cv_type.eq.'circular') then
 cvv(curve_temp)%begin%next%previous=>cvv(curve_temp)%eend%previous
 cvv(curve_temp)%eend%previous%next=>cvv(curve_temp)%begin%next
end if

! call check_list_c(1,'down  ')

call deletee(cvv(curve_number))
cvv(curve_number)=cvv(curve_temp)

cvv(curve_temp)%status='asleep'
cvv(curve_temp)%cv_type='   '
cvv(curve_temp)%begin_end=0
cvv(curve_temp)%end_end=0
cvv(curve_temp)%total=0
cvv(curve_temp)%wave=0
nullify(cvv(curve_temp)%begin); nullify(cvv(curve_temp)%eend)


contains


subroutine insertt

implicit none

!type(geo_info) :: g_cell, g_new
!type(adss_info) :: ads, ads_new

nullify(temp_tp%previous); nullify(temp_tp%next)

if(associated(current%previous)) then
! The insert point is in the middle of the discontinuity curve.   
 temp_tp%previous=>current%previous
 temp_tp%next=>current
 current%previous%next=>temp_tp
 current%previous=>temp_tp
else
 if(associated(cvv(curve_temp)%begin%next)) then
! The insert point is in front of the first critical cell in 
! the curve. 
  temp_tp%next=>current
  current%previous=>temp_tp
  cvv(curve_temp)%begin%next=>temp_tp
 else
! The curve is still empty and 'temp_new' is going to be the 
! first member of the curve. 'temp' is the curve's 'begin'.
  cvv(curve_temp)%begin%next=>temp_tp
  cvv(curve_temp)%eend%previous=>temp_tp
 end if
end if

end subroutine insertt

 
end subroutine revers_curve


end module revers_discontinuity_curves