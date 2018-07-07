module smoothen_curves
! The ordinary states of some discontinuity cells may exceed the left or 
! right states at the end of a step of computation. If this happens, 
! the exceeded parts of the ordinary states will be cut of and given to 
! the neighboring critical cells. This module is for this purpose.

! The current one is only for interfaces.

use discontinuity_curves

implicit none

public smoothen_cv


contains


subroutine smoothen_cv

implicit none

call truncation

call distribution


contains


subroutine truncation

implicit none

type(critical_cell), pointer :: tempx
type(geo_info) :: g_cell
type(phy_info) :: p_cell
real(8) :: density_difference
integer :: head_mark, i
logical :: head_switch

type(state) :: l_state, r_state
real(8) :: density_l, density_r, density_or
real(8) :: x_velo, y_velo, x_velo_p, y_velo_p, x_velo_n, y_velo_n
real(8), dimension(2) :: normal

do i=1, curves_number
 if(cvv(i)%status.eq.'awake ') then
  tempx=>cvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
   
  do while(associated(tempx).and.head_switch) 
   g_cell=tempx%g_cell; p_cell=tempx%p_cell
   tempx%l_store=0.0d0; tempx%r_store=0.0d0

   normal=(/1.0d0,0.0d0/)
   call find_physical_state(p_cell%l_state,normal,l_state)
   call find_physical_state(p_cell%r_state,normal,r_state)

   density_or=density(p_cell%or_state)
   density_l=l_state%value(1); density_r=r_state%value(1)

   x_velo=0.5d0*(l_state%value(2)+r_state%value(2))
   y_velo=0.5d0*(l_state%value(3)+r_state%value(3))
   x_velo_p=0.5d0*(x_velocity(tempx%previous%p_cell%l_state)+ &
            x_velocity(tempx%previous%p_cell%r_state))
   y_velo_p=0.5d0*(y_velocity(tempx%previous%p_cell%l_state)+ &
            y_velocity(tempx%previous%p_cell%r_state))
   x_velo_n=0.5d0*(x_velocity(tempx%next%p_cell%l_state)+ &
            x_velocity(tempx%next%p_cell%r_state))
   y_velo_n=0.5d0*(y_velocity(tempx%next%p_cell%l_state)+ &
            y_velocity(tempx%next%p_cell%r_state))

   density_difference=0.0d0

   if(density_l.lt.density_r) then
    if(density_or.lt.density_l) then
     density_difference=density_l-density_or
    else
     if(density_or.gt.density_r) then
      density_difference=density_r-density_or
     end if
    end if
   else
    if(density_or.gt.density_l) then
     density_difference=density_l-density_or
    else
     if(density_or.lt.density_r) then
      density_difference=density_r-density_or
     end if
    end if
   end if
     
   if(density_difference.gt.0.0d0) then
     
    tempx%l_store%value(1)=density_difference
	tempx%l_store%value(2)=density_difference*x_velo_p
	tempx%l_store%value(3)=density_difference*y_velo_p
    tempx%l_store%value(4)=0.5d0*density_difference*(x_velo_p*x_velo_p+y_velo_p*y_velo_p)
     
    tempx%r_store%value(1)=density_difference
	tempx%r_store%value(2)=density_difference*x_velo_n
	tempx%r_store%value(3)=density_difference*y_velo_n
    tempx%r_store%value(4)=0.5d0*density_difference*(x_velo_n*x_velo_n+y_velo_n*y_velo_n)
      
   else
     
    if(density_difference.lt.0.0d0) then
     
     tempx%l_store%value(1)=density_difference
	 tempx%l_store%value(2)=density_difference*x_velo
	 tempx%l_store%value(3)=density_difference*y_velo
     tempx%l_store%value(4)=0.5d0*density_difference*(x_velo*x_velo+y_velo*y_velo)
     
     tempx%r_store%value(1)=density_difference
	 tempx%r_store%value(2)=density_difference*x_velo
	 tempx%r_store%value(3)=density_difference*y_velo
     tempx%r_store%value(4)=0.5d0*density_difference*(x_velo*x_velo+y_velo*y_velo)   
    
	end if
   end if
   
   tempx=>tempx%next
   if(associated(tempx)) then
    head_switch=(tempx%address%idx.ne.head_mark)
   else
    exit
   end if
  end do
  
 end if
end do

end subroutine truncation


subroutine distribution

implicit none

type(critical_cell), pointer :: tempx
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(state) :: p_difference, n_difference
integer :: head_mark, i
logical :: head_switch

do i=1, curves_number
 if(cvv(i)%status.eq.'awake ') then
  tempx=>cvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
  
  do while(associated(tempx).and.head_switch) 
   g_cell=tempx%g_cell; p_cell=tempx%p_cell
   p_difference=tempx%l_store; n_difference=tempx%r_store
   if(associated(tempx%previous)) then
    if(associated(tempx%next)) then
     tempx%previous%p_cell%or_state=tempx%previous%p_cell%or_state-0.5d0*p_difference
     tempx%next%p_cell%or_state=tempx%next%p_cell%or_state-0.5d0*n_difference
     tempx%p_cell%or_state=tempx%p_cell%or_state+0.5d0*(p_difference+n_difference)
    else
     tempx%previous%p_cell%or_state=tempx%previous%p_cell%or_state-0.5d0*p_difference
     tempx%p_cell%or_state=tempx%p_cell%or_state+0.5d0*p_difference
    end if
   else
    if(associated(tempx%next)) then
     tempx%next%p_cell%or_state=tempx%next%p_cell%or_state-0.5d0*n_difference
     tempx%p_cell%or_state=tempx%p_cell%or_state+0.5d0*n_difference
    end if
   end if						 	 									 								     
   
   tempx=>tempx%next
   if(associated(tempx)) then
    head_switch=(tempx%address%idx.ne.head_mark)
   else
    exit
   end if
  end do
  
 end if
end do
   
end subroutine distribution


end subroutine smoothen_cv


end module smoothen_curves