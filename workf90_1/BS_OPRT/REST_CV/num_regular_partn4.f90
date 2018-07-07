module numerical_regular_partners
! This module is for the numerical regularization on discontinuity curves.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use geo_info_cell  
! gif_cell.f90'

implicit none

public  numerical_regularize_p, make_partner_memory
private clear_numerical_regular_p, compute_states, distribute_states


contains


subroutine make_partner_memory

implicit none

integer :: i

do i=1,cvn
 if(acvv(i)%status.eq.'awake ') then
  call partner_memory(acvv(i))
 end if
end do  


contains


subroutine partner_memory(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end

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
 address=temp%address

 temp%partner_memory=temp%partner   

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine partner_memory


end subroutine make_partner_memory


subroutine numerical_regularize_p
! This subroutine performs the numerical regularization.

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call clear_numerical_regular_p(acvv(i))
 end if
end do

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call compute_states(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call distribute_states(acvv(i))
 end if
end do

end subroutine numerical_regularize_p


subroutine clear_numerical_regular_p(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end

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
 address=temp%address

 temp%f_cell%numerical_regularization(1)=0.0d0   

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine clear_numerical_regular_p


subroutine compute_states(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell, g_cellp
type(phy_info) :: p_cell, p_cellp

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

real(8) :: l_area, r_area, length, l_areap, r_areap, lengthp
real(8) :: uu_l, vv_l, uu_r, vv_r, uu_or_temp, vv_or_temp, kin_ee
real(8) :: uu_lp, vv_lp, uu_rp, vv_rp, uu_or_tempp, vv_or_tempp, kin_eep
type(state) :: or_auxi, or_state_temp, or_state_tempp, or_state_auxi, or_state_auxip, or_total
type(state) :: physical_difference, physical_differencep
real(8) :: alpha, beta

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
 address=temp%address
 p_cell=temp%p_cell; g_cell=temp%g_cell
   
 if(temp%partner_memory.eq.'next    '.and.temp%next%partner_memory.eq.'previous') then

  g_cellp=temp%next%g_cell; p_cellp=temp%next%p_cell
   
  if(g_cell%g_type.eq.'xy '.or.g_cellp%g_type.eq.'xy ')  then
   
   l_area=side_area(g_cell,'left  '); r_area=1.0d0-l_area
   length=length_of_curve_segment(g_cell)
   l_areap=side_area(g_cellp,'left  '); r_areap=1.0d0-l_areap 
   lengthp=length_of_curve_segment(g_cellp)   
   
   or_total=p_cell%or_state+p_cellp%or_state
   or_state_temp=l_area*p_cell%l_state+r_area*p_cell%r_state
   or_state_tempp=l_areap*p_cellp%l_state+r_areap*p_cellp%r_state
   
   alpha=lengthp/(length+lengthp); beta=1.0d0-alpha
   or_state_auxi=alpha*or_state_temp+beta*(or_total-or_state_tempp)
   or_state_auxip=alpha*(or_total-or_state_temp)+beta*or_state_tempp  
   
   uu_l=x_velocity(p_cell%l_state); uu_r=x_velocity(p_cell%r_state)
   vv_l=y_velocity(p_cell%l_state); vv_r=y_velocity(p_cell%r_state)
   uu_or_temp=l_area*uu_l+r_area*uu_r; vv_or_temp=l_area*vv_l+r_area*vv_r
    
   or_state_temp%value(1)=or_state_auxi%value(1)
   or_state_temp%value(2)=or_state_auxi%value(1)*uu_or_temp
   or_state_temp%value(3)=or_state_auxi%value(1)*vv_or_temp
   or_state_temp%value(4)=0.5d0*or_state_auxi%value(1)*(uu_or_temp**2.0d0+vv_or_temp**2.0d0)
   
   kin_ee=or_state_auxi%value(4)-internal_energy(or_state_auxi)   
   
   uu_lp=x_velocity(p_cellp%l_state); uu_rp=x_velocity(p_cellp%r_state)
   vv_lp=y_velocity(p_cellp%l_state); vv_rp=y_velocity(p_cellp%r_state)
   uu_or_tempp=l_areap*uu_lp+r_areap*uu_rp; vv_or_tempp=l_areap*vv_lp+r_areap*vv_rp
   
   or_state_tempp%value(1)=or_state_auxip%value(1)
   or_state_tempp%value(2)=or_state_auxip%value(1)*uu_or_tempp
   or_state_tempp%value(3)=or_state_auxip%value(1)*vv_or_tempp
   or_state_tempp%value(4)=0.5d0*or_state_auxip%value(1)*(uu_or_tempp**2.0d0+vv_or_tempp**2.0d0)
   
   kin_eep=or_state_auxip%value(4)-internal_energy(or_state_auxip)      
   
   physical_difference=or_state_auxi-or_state_temp
   physical_difference%value(4)=kin_ee-or_state_temp%value(4)
   physical_differencep=or_state_auxip-or_state_tempp    
   physical_differencep%value(4)=kin_eep-or_state_tempp%value(4)
   
   or_state_auxi=or_state_auxi-alpha*physical_difference+beta*physical_differencep
   or_state_auxip=or_state_auxip+alpha*physical_difference-beta*physical_differencep 

   or_state_auxi%gamma=error_data; or_state_auxi%gamma=error_data
     
   temp%f_cell%numerical_regularization(1)=or_state_auxi
   temp%next%f_cell%numerical_regularization(1)=or_state_auxip
    
!   temp=>temp%next
    
  end if
     
 end if
   
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine compute_states


subroutine distribute_states(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell, g_cellp
type(phy_info) :: p_cell, p_cellp
type(flux_info) :: f_cell, f_cellp

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
 address=temp%address
 p_cell=temp%p_cell; g_cell=temp%g_cell; f_cell=temp%f_cell
  
 if(temp%partner_memory.eq.'next    '.and.temp%next%partner_memory.eq.'previous') then
  p_cellp=temp%next%p_cell; g_cellp=temp%next%g_cell; f_cellp=temp%next%f_cell
  if(g_cell%g_type.eq.'xy '.or.g_cellp%g_type.eq.'xy ') then
   p_cell%or_state=f_cell%numerical_regularization(1)
   p_cellp%or_state=f_cellp%numerical_regularization(1)
   temp%next%p_cell=p_cellp
   temp%p_cell=p_cell
  end if
!  temp=>temp%next
 end if
       
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine distribute_states


end module numerical_regular_partners