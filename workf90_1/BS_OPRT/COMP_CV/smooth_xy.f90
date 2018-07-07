module smoothen_xy_cell_averages
! After the computation on discontinuity curves, the ordinary states of some 'xy'-type 
! critical cells may exceed the ranges between the left and right states. If this happens, 
! the exceeded parts of the ordinary states will be cut and given to the neighboring 
! critical cells. This module is for this purpose.

use auxiliary_discontinuity_curves

implicit none

public smoothen_xy


contains


subroutine smoothen_xy

implicit none

type(auxi_crit_cell), pointer :: tempx
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(state) :: difference
integer :: head_mark, i, single
logical :: head_switch

do i=1, curves_number
 if(acvv(i)%status.eq.'awake ') then
  tempx=>acvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
  do while(associated(tempx).and.head_switch) 
   g_cell=tempx%g_cell; p_cell=tempx%p_cell
   difference=0.d0
   if(g_cell%g_type.eq.'xy ') then
    call find_single(g_cell,single)
    select case(g_cell%point(single))
     case('left  ')
      call side_cut(p_cell%l_state,p_cell%r_state,p_cell%or_state, &
                    difference)
     case('right ')
      call side_cut(p_cell%r_state,p_cell%l_state,p_cell%or_state, &
                    difference)
     case default; call error_message; pause
    end select 
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


end subroutine smoothen_xy

end module smoothen_xy_cell_averages