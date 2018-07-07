module packed_separation
! interface segments in packed discontinuity cells may go across each other. This model is to keep them
! spearated by revising the ordinary states of the two packed discontinuity cells. 

! This treatment is provisional. Better treatment should be to make topology change of the interface.

use discontinuity_curves

implicit none

public packed_separate


contains


subroutine packed_separate

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

type(critical_cell), pointer :: tempxn
real(8) :: density_l, density_r, density_or, x_velo, y_velo, pressure
real(8) :: area_this, area_that, area_exceeded

do i=1, curves_number
 if(cvv(i)%status.eq.'awake ') then
  tempx=>cvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
   
  do while(associated(tempx).and.head_switch) 
    
   if(tempx%l_stk.ne.adss_info(0,0).or.tempx%r_stk.ne.adss_info(0,0)) then
    g_cell=tempx%g_cell; p_cell=tempx%p_cell
     
!    if(tempx%l_stk.ne.adss_info(0,0).and.tempx%r_stk.ne.adss_info(0,0)) call error_message
    
    density_l=density(p_cell%l_state) 
    density_r=density(p_cell%r_state)   
    density_or=density(p_cell%or_state)
    
    if(tempx%l_stk.ne.adss_info(0,0)) then
      
     area_this=(density_or-density_l)/(density_r-density_l)
	 call visit(tempx%l_stk,tempxn)
     if(tempxn%l_stk.eq.tempx%address) then
	  if(dabs(density(tempxn%p_cell%l_state)-density_l).gt.1.0d-12) call error_message
      area_that=density(tempxn%p_cell%or_state)-density(tempxn%p_cell%l_state)
      area_that=area_that/(density(tempxn%p_cell%r_state)-density(tempxn%p_cell%l_state))
     else
	  if(dabs(density(tempxn%p_cell%r_state)-density_l).gt.1.0d-12) call error_message
      area_that=density(tempxn%p_cell%or_state)-density(tempxn%p_cell%r_state)
      area_that=area_that/(density(tempxn%p_cell%l_state)-density(tempxn%p_cell%r_state))      
     end if
     
    else
     
     area_this=(density_or-density_r)/(density_l-density_r)
	 call visit(tempx%r_stk,tempxn)
     if(tempxn%l_stk.eq.tempx%address) then
	  if(dabs(density(tempxn%p_cell%l_state)-density_r).gt.1.0d-12) call error_message
      area_that=density(tempxn%p_cell%or_state)-density(tempxn%p_cell%l_state)
      area_that=area_that/(density(tempxn%p_cell%r_state)-density(tempxn%p_cell%l_state))
     else
	  if(dabs(density(tempxn%p_cell%r_state)-density_r).gt.1.0d-12) call error_message
      area_that=density(tempxn%p_cell%or_state)-density(tempxn%p_cell%r_state)
      area_that=area_that/(density(tempxn%p_cell%l_state)-density(tempxn%p_cell%r_state))      
     end if
     
	end if
     
    area_exceeded=0.5d0*dmax1(area_this+area_that-1.0d0+1.0d-1,0.0d0)
    x_velo=0.5d0*(x_velocity(p_cell%l_state)+x_velocity(p_cell%r_state))
    y_velo=0.5d0*(x_velocity(p_cell%l_state)+y_velocity(p_cell%r_state))
          
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

     difference=0.0d0
     
     if(tempx%l_stk.ne.adss_info(0,0)) then
	  difference%value(1)=area_exceeded*(density_r-density_l)
     else
	  difference%value(1)=area_exceeded*(density_l-density_r)
     end if
     difference%value(2)=difference%value(1)*x_velo
     difference%value(3)=difference%value(1)*y_velo
     difference%value(4)=0.5d0*difference%value(1)*(x_velo*x_velo+y_velo*y_velo)  
    
! This is only for $\gamma$-law gases.
    
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   
        
    tempx%memo=difference
    
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
type(state) :: difference
integer :: head_mark, i
logical :: head_switch

do i=1, curves_number
 if(cvv(i)%status.eq.'awake ') then
  tempx=>cvv(i)%begin%next
  head_mark=tempx%address%idx; head_switch=.true.
  
  do while(associated(tempx).and.head_switch) 

   if(tempx%l_stk.ne.adss_info(0,0).or.tempx%r_stk.ne.adss_info(0,0)) then
    g_cell=tempx%g_cell; p_cell=tempx%p_cell
    difference=tempx%memo; tempx%memo=0.0d0
    tempx%p_cell%or_state=tempx%p_cell%or_state-difference
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


end subroutine packed_separate


end module packed_separation