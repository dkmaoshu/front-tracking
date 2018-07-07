module numerical_regularization
! This module is for the numerical regularization on discontinuity curves.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use geo_info_cell  
! gif_cell.f90'

implicit none

public  numerical_regularize
private clear_numerical_regularization, compute_den_inte_differences, compute_mom_kin_differences, & 
        distribute_differences1, distribute_differences2


contains


subroutine numerical_regularize
! This subroutine performs the numerical regularization.

implicit none

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call clear_numerical_regularization(acvv(i))
 end if
end do

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call compute_den_inte_differences(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call distribute_differences1(acvv(i))
 end if
end do

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call compute_mom_kin_differences(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call distribute_differences2(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

end subroutine numerical_regularize


subroutine clear_numerical_regularization(acv)

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

 temp%f_cell%numerical_regularization(1)=error_data   
 temp%f_cell%numerical_regularization(2)=error_data

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine clear_numerical_regularization


subroutine compute_den_inte_differences(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

real(8) :: l_area, r_area, or_den_temp, or_inte_temp, l_inte, r_inte
type(state) :: physical_difference
integer :: kkk

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
    
 select case(g_cell%g_type)
    
  case('xx ','yy')
   select case(g_cell%point(1))
    case('left  ')
	 l_area=0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2))+0.5d0
	 r_area=1.0d0-l_area
    case('right ')
	 l_area=0.5d0-0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2))
	 r_area=1.0d0-l_area
    case default; call error_message
   end select
         
  case('xy ')
   l_area=side_area(g_cell,'left  ')
   r_area=1.0d0-l_area

  case default; call error_message
   
 end select

 or_den_temp=l_area*p_cell%l_state%value(1)+r_area*p_cell%r_state%value(1)
 l_inte=internal_energy(p_cell%l_state); r_inte=internal_energy(p_cell%r_state)
 or_inte_temp=l_area*l_inte+r_area*r_inte

 physical_difference%value(1)=p_cell%or_state%value(1)-or_den_temp
 do kkk=2, 3
  physical_difference%value(kkk)=0.0d0
 end do
 physical_difference%value(4)=internal_energy(p_cell%or_state)-or_inte_temp


 select case(temp%partner)

  case('single  ')
   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   
  case('previous')
   temp%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   temp%previous%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   
  case('next    ')
   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%next%f_cell%numerical_regularization(1)=0.5d0*physical_difference

  case default; call error_message

 end select

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine compute_den_inte_differences


subroutine compute_mom_kin_differences(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

real(8) :: uu_l, uu_r, vv_l, vv_r,kin_ee
real(8) :: den_or_temp, uu_or_temp, vv_or_temp
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
   
 uu_l=x_velocity(p_cell%l_state); uu_r=x_velocity(p_cell%r_state)
 vv_l=y_velocity(p_cell%l_state); vv_r=y_velocity(p_cell%r_state)
 kin_ee=p_cell%or_state%value(4)-internal_energy(p_cell%or_state)

 den_or_temp=p_cell%or_state%value(1)
 uu_or_temp=0.5d0*(uu_l+uu_r)
 vv_or_temp=0.5d0*(vv_l+vv_r)

 or_state_temp%value(1)=den_or_temp
 or_state_temp%value(2)=den_or_temp*uu_or_temp
 or_state_temp%value(3)=den_or_temp*vv_or_temp
 or_state_temp%value(4)=0.5d0*den_or_temp*(uu_or_temp**2.0d0+vv_or_temp**2.0d0)

 physical_difference=p_cell%or_state-or_state_temp
 physical_difference%value(4)=kin_ee-or_state_temp%value(4)
     
 select case(temp%partner)

  case('single  ')
   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   
  case('previous')
   temp%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   temp%previous%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   
  case('next    ')
   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%next%f_cell%numerical_regularization(1)=0.5d0*physical_difference

  case default; call error_message

 end select
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine compute_mom_kin_differences


subroutine distribute_differences1(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(flux_info) :: f_cell, f_cellp, f_celln

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
 f_cellp=temp%previous%f_cell; f_celln=temp%next%f_cell

 select case(temp%partner)

  case('single  ')
   p_cell%or_state=p_cell%or_state-2.0d0*rr*f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-2.0d0*rr*f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+2.0d0*rr*temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+2.0d0*rr*temp%next%f_cell%numerical_regularization(1)

  case('previous')
   p_cell%or_state=p_cell%or_state-2.0d0*rr*temp%previous%f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-2.0d0*rr*f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+2.0d0*rr*temp%previous%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+2.0d0*rr*temp%next%f_cell%numerical_regularization(1)

  case('next    ')
   p_cell%or_state=p_cell%or_state-2.0d0*rr*f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-2.0d0*rr*temp%next%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+2.0d0*rr*temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+2.0d0*rr*temp%next%next%f_cell%numerical_regularization(1)

  case default; call error_message

 end select

 temp%p_cell=p_cell
       
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine distribute_differences1


subroutine distribute_differences2(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(flux_info) :: f_cell, f_cellp, f_celln

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
 f_cellp=temp%previous%f_cell; f_celln=temp%next%f_cell

 select case(temp%partner)

  case('single  ')
   p_cell%or_state=p_cell%or_state-f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+temp%next%f_cell%numerical_regularization(1)

  case('previous')
   p_cell%or_state=p_cell%or_state-temp%previous%f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+temp%previous%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+temp%next%f_cell%numerical_regularization(1)

  case('next    ')
   p_cell%or_state=p_cell%or_state-f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-temp%next%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+temp%next%next%f_cell%numerical_regularization(1)

  case default; call error_message

 end select

 temp%p_cell=p_cell
       
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine distribute_differences2


end module numerical_regularization