module numerical_regular_partners
! This module is for the numerical regularization on discontinuity curves.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use geo_info_cell  
! gif_cell.f90'

implicit none

public  numerical_regularize_p, make_partner_memory
private clear_numerical_regular_p, compute_differences_p, distribute_differences_p


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
  call compute_differences_p(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call distribute_differences_p(acvv(i))
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


subroutine compute_differences_p(acv)


type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

real(8) :: l_area, r_area, uu_l, uu_r, vv_l, vv_r, ee_l, ee_r
real(8) :: den_or_temp, uu_or_temp, vv_or_temp, ee_or_temp
type(state) :: or_state_temp, physical_difference

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

 if(temp%partner_memory.eq.'previous'.or.temp%partner_memory.eq.'next    ') then
    
  if(g_cell%g_type.eq.'xy ')  then
    
   l_area=side_area(g_cell,'left  ')
   r_area=1.0d0-l_area

   uu_l=x_velocity(p_cell%l_state); uu_r=x_velocity(p_cell%r_state)
   vv_l=y_velocity(p_cell%l_state); vv_r=y_velocity(p_cell%r_state)
   ee_l=internal_energy(p_cell%l_state); ee_r=internal_energy(p_cell%r_state)

   den_or_temp=p_cell%or_state%value(1)
   uu_or_temp=l_area*uu_l+r_area*uu_r
   vv_or_temp=l_area*vv_l+r_area*vv_r
   ee_or_temp=l_area*ee_l+r_area*ee_r

   or_state_temp%value(1)=den_or_temp
   or_state_temp%value(2)=den_or_temp*uu_or_temp
   or_state_temp%value(3)=den_or_temp*vv_or_temp
   or_state_temp%value(4)=0.5d0*den_or_temp*(uu_or_temp**2.0d0+vv_or_temp**2.0d0)+ee_or_temp

   physical_difference=p_cell%or_state-or_state_temp
    
!   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%f_cell%numerical_regularization(1)=physical_difference
    
  end if
  
 end if
   
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine compute_differences_p


subroutine distribute_differences_p(acv)


type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(flux_info) :: f_cell

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
  
 if(g_cell%g_type.eq.'xy ') then
  
  select case(temp%partner_memory)
   
   case('previous')
    p_cell%or_state=p_cell%or_state-f_cell%numerical_regularization(1)
    p_cell%or_state=p_cell%or_state+temp%previous%f_cell%numerical_regularization(1)

   case('next    ')
    p_cell%or_state=p_cell%or_state-f_cell%numerical_regularization(1)
    p_cell%or_state=p_cell%or_state+temp%next%f_cell%numerical_regularization(1)

  end select

 end if

 temp%p_cell=p_cell
       
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine distribute_differences_p


end module numerical_regular_partners