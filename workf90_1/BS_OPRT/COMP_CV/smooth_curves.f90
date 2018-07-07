module smoothen_curves
! The ordinary states of some discontinuity cells may exceed the left or 
! right states at the end of a step of computation. If this happens, 
! the exceeded parts of the ordinary states will be cut of and given to 
! the neighboring critical cells. This module is for this purpose.

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
type(state) :: difference
integer :: head_mark, i
logical :: head_switch

do i=1, curves_number
 if(cvv(i)%status.eq.'awake ') then
  tempx=>cvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
   
  do while(associated(tempx).and.head_switch) 
   g_cell=tempx%g_cell; p_cell=tempx%p_cell
   call side_cut(p_cell%l_state,p_cell%r_state,p_cell%or_state, &
                 p_cell%wv_nb,difference)

!    if(dabs(difference%value(1)).gt.1.0d-12) then
!	 print*, 'I=', g_cell%x_idx, 'J=', g_cell%y_idx, difference%value(1)
!    end if

   tempx%memo=difference
   call side_cut(p_cell%r_state,p_cell%l_state,p_cell%or_state, &
                 p_cell%wv_nb,difference)

!    if(dabs(difference%value(1)).gt.1.0d-12) then
!	 print*, 'I=', g_cell%x_idx, 'J=', g_cell%y_idx, difference%value(1)
!    end if

   tempx%memo=tempx%memo+difference 
  
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
type(state) :: difference
integer :: head_mark, i
logical :: head_switch

do i=1, curves_number
 if(cvv(i)%status.eq.'awake ') then
  tempx=>cvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
  
  do while(associated(tempx).and.head_switch) 
   g_cell=tempx%g_cell; p_cell=tempx%p_cell
   difference=tempx%memo; tempx%memo=0.0d0
   if(associated(tempx%previous)) then
    if(associated(tempx%next)) then
     tempx%previous%p_cell%or_state=tempx%previous%p_cell%or_state+ &
                                    0.5d0*difference
     tempx%next%p_cell%or_state=tempx%next%p_cell%or_state+  &
                                0.5d0*difference
     tempx%p_cell%or_state=tempx%p_cell%or_state-difference
    else
     tempx%previous%p_cell%or_state=tempx%previous%p_cell%or_state+ &
                                    0.5d0*difference
     tempx%p_cell%or_state=tempx%p_cell%or_state-0.5d0*difference
    end if
   else
    if(associated(tempx%next)) then
     tempx%next%p_cell%or_state=tempx%next%p_cell%or_state+ &
                                0.5d0*difference
     tempx%p_cell%or_state=tempx%p_cell%or_state-0.5d0*difference
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