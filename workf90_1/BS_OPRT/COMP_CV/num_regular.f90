module numerical_regularization
! This module is for the numerical regularization on discontinuity curves.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use geo_info_cell  
! gif_cell.f90'

implicit none

public  numerical_regularize
private clear_numerical_regularization, compute_differences, distribute_differences


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
  call compute_differences(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call distribute_differences(acvv(i))
 end if
end do

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


subroutine compute_differences(acv)


type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

real(8) :: area_difference, single_area
type(state) :: physical_difference

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
    
  case('xx ')
   area_difference=g_cell%mdis_posi(1)-0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2))
   select case(g_cell%point(1))
    case('left  ')
	 physical_difference=area_difference*(p_cell%l_state-p_cell%r_state)
    case('right ')
	 physical_difference=area_difference*(p_cell%r_state-p_cell%l_state)
    case default; call error_message
   end select
         
  case('yy ')
   area_difference=g_cell%mdis_posi(2)-0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2))
   select case(g_cell%point(1)) 
    case('left  ')
	 physical_difference=area_difference*(p_cell%l_state-p_cell%r_state)
    case('right ')
	 physical_difference=area_difference*(p_cell%r_state-p_cell%l_state)
    case default; call error_message
   end select
    
  case('xy ')
   single_area=side_area(g_cell,'left  ')
   area_difference=g_cell%mdis_posi(1)-single_area
   physical_difference=area_difference*(p_cell%l_state-p_cell%r_state)

  case default; call error_message
   
 end select

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

end subroutine compute_differences


subroutine distribute_differences(acv)


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

 select case(temp%partner)

  case('single  ')
   p_cell%or_state%value(1)=p_cell%or_state%value(1)-0.5d0*r*f_cell%numerical_regularization(1)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)-0.5d0*r*f_cell%numerical_regularization(2)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)+0.5d0*r*temp%previous%f_cell%numerical_regularization(2)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)+0.5d0*r*temp%next%f_cell%numerical_regularization(1)%value(1)

  case('previous')
   p_cell%or_state%value(1)=p_cell%or_state%value(1)-0.5d0*r*temp%previous%f_cell%numerical_regularization(1)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)-0.5d0*r*f_cell%numerical_regularization(2)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)+0.5d0*r*temp%previous%previous%f_cell%numerical_regularization(2)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)+0.5d0*r*temp%next%f_cell%numerical_regularization(1)%value(1)

  case('next    ')
   p_cell%or_state%value(1)=p_cell%or_state%value(1)-0.5d0*r*f_cell%numerical_regularization(1)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)-0.5d0*r*temp%next%f_cell%numerical_regularization(2)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)+0.5d0*r*temp%previous%f_cell%numerical_regularization(2)%value(1)
   p_cell%or_state%value(1)=p_cell%or_state%value(1)+0.5d0*r*temp%next%next%f_cell%numerical_regularization(1)%value(1)

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

end subroutine distribute_differences


end module numerical_regularization