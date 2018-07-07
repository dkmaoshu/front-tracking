module numerical_regularization
! This module is for the numerical regularization on discontinuity curves.

! THIS IS THE SIX TH VERSION OF THE FILE.

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

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call clear_numerical_regularization(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

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

! call check_list_ac(1,'down  ')

end subroutine numerical_regularize


subroutine clear_numerical_regularization(acv)

type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

nullify(temp)
if(associated(acv%begin_boundary%next)) then
 temp=>acv%begin_boundary%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend_boundary%previous%address%idx
 end if
end if

do while(associated(temp).and.head_switch) 
 address=temp%address

 temp%f_cell%numerical_regularization(1)=error_data   
 temp%f_cell%numerical_regularization(2)=error_data
 temp%f_cell%numer_reg_coefficient=0.0d0

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
type(geo_info) :: g_cell, g_cell_p
type(phy_info) :: p_cell, p_cell_p

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

real(8) :: l_area, r_area, l_area_p, r_area_p, length, length_p, num_reg_coef, mdposi, conser_area
type(state) :: physical_difference, or_state_temp, ttt

nullify(temp)
if(associated(acv%begin_boundary%next)) then
 temp=>acv%begin_boundary%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend_boundary%previous%address%idx
 end if
end if
   
do while(associated(temp).and.head_switch) 
 address=temp%address
 p_cell=temp%p_cell; g_cell=temp%g_cell; num_reg_coef=0.0d0
   
 select case(temp%partner)
   
  case('single  ')
    
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
   
   or_state_temp=l_area*p_cell%l_state+r_area*p_cell%r_state
   physical_difference=p_cell%or_state-or_state_temp
    
! Modification in the sixth verion over the fifth version.
!##########################################################################
   select case(g_cell%g_type)
    case('xx ') 
     select case(g_cell%point(1))
      case('left  '); num_reg_coef=g_cell%mdis_posi(1)-l_area+0.5d0    
      case('right '); num_reg_coef=g_cell%mdis_posi(1)-r_area+0.5d0
      case default; call error_message
     end select
    case('yy ') 
     select case(g_cell%point(1))
      case('left  '); num_reg_coef=g_cell%mdis_posi(2)-l_area+0.5d0
      case('right '); num_reg_coef=g_cell%mdis_posi(2)-r_area+0.5d0
     end select
    case('xy ')
     conser_area=(p_cell%l_state%value(1)-p_cell%r_state%value(1))
	 conser_area=(p_cell%or_state%value(1)-p_cell%r_state%value(1))/conser_area
     num_reg_coef=conser_area-l_area
    case default; call error_message
   end select
!##########################################################################
   
  case('previous', 'next    ')
    
   if(temp%partner.eq.'previous') then
    g_cell_p=temp%previous%g_cell; p_cell_p=temp%previous%p_cell
   else
    g_cell_p=temp%next%g_cell; p_cell_p=temp%next%p_cell
   end if	 
     
   call side_area_length_partnered(g_cell,g_cell_p,l_area,r_area,length)
   call side_area_length_partnered(g_cell_p,g_cell,l_area_p,r_area_p,length_p)
   or_state_temp=l_area*p_cell%l_state+r_area*p_cell%r_state
   or_state_temp=or_state_temp+l_area_p*p_cell_p%l_state+r_area_p*p_cell_p%r_state
   
   select case(g_cell_p%x_idx+g_cell_p%y_idx-g_cell%x_idx-g_cell%y_idx) 
    case(1)
     if(g_cell%point(1).eq.'left  ') then
      ttt=p_cell%or_state+p_cell_p%r_state
     else
      ttt=p_cell%or_state+p_cell_p%l_state
     end if
    case(-1)
     if(g_cell%point(3).eq.'left  ') then
      ttt=p_cell%or_state+p_cell_p%r_state
     else
      ttt=p_cell%or_state+p_cell_p%l_state
     end if
    case default; call error_message
   end select		 		 	  	 	  	   
    
   physical_difference=ttt-or_state_temp
     
! Modification in the sixth verion over the fifth version.
!##########################################################################
   select case(g_cell_p%x_idx+g_cell_p%y_idx-g_cell%x_idx-g_cell%y_idx) 
    case(1)
 !    select case(g_cell%point(1))
 !     case('left  ') 
       select case(g_cell%g_type)
        case('xx '); mdposi=g_cell%mdis_posi(1)
        case('yy '); mdposi=g_cell%mdis_posi(2)
        case default; call error_message
       end select
 !     case('right ') 
 !      select case(g_cell%g_type)
 !       case('xx '); mdposi=g_cell%mdis_posi(1)
 !       case('yy '); mdposi=g_cell%mdis_posi(2)
 !       case default; call error_message
 !      end select
 !     case default; call error_message
 !    end select
    case(-1)
 !    select case(g_cell%point(1))
 !     case('left  ') 
       select case(g_cell%g_type)
        case('xx '); mdposi=g_cell%mdis_posi(1)+1.0d0
        case('yy '); mdposi=g_cell%mdis_posi(2)+1.0d0
        case default; call error_message
       end select
 !     case('right ') 
 !      select case(g_cell%g_type)
 !       case('xx '); mdposi=g_cell%mdis_posi(1)+1.0d0
 !       case('yy '); mdposi=g_cell%mdis_posi(2)+1.0d0
 !       case default; call error_message
 !      end select
 !     case default; call error_message
 !    end select
    case default; call error_message
   end select 
     
   select case(g_cell%point(1))
    case('left  '); num_reg_coef=mdposi-l_area-l_area_p+0.5d0
    case('right '); num_reg_coef=mdposi-r_area-r_area_p+0.5d0
    case default; call error_message
   end select
!##########################################################################
   
  case default; call error_message
   
 end select 
        
 select case(temp%partner)
     
  case('single  ')
   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   temp%f_cell%numer_reg_coefficient=num_reg_coef
   
  case('previous')
   temp%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   temp%f_cell%numer_reg_coefficient=num_reg_coef
   temp%previous%f_cell%numerical_regularization(2)=0.5d0*physical_difference
   temp%previous%f_cell%numer_reg_coefficient=num_reg_coef
   
  case('next    ')
   temp%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%f_cell%numer_reg_coefficient=num_reg_coef
   temp%next%f_cell%numerical_regularization(1)=0.5d0*physical_difference
   temp%next%f_cell%numer_reg_coefficient=num_reg_coef
    
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
type(flux_info) :: f_cell, f_cellp, f_celln

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch
real(8) :: coeff_current, coeff_previous, coeff_next, coeff_previous_previous, coeff_next_next

!coeff=0.0d0
!coeff=rr
!coeff=0.5d0*rr

!coeff=0.01d0*rr

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

 coeff_current=dabs(f_cell%numer_reg_coefficient)
 coeff_previous=dabs(temp%previous%f_cell%numer_reg_coefficient)
 coeff_next=dabs(temp%next%f_cell%numer_reg_coefficient)
 coeff_previous_previous=0.0d0; coeff_next_next=0.0d0
 
 select case(temp%partner)

  case('single  ')
   p_cell%or_state=p_cell%or_state-0.5d0*coeff_current*f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-0.5d0*coeff_current*f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+0.5d0*coeff_previous*temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+0.5d0*coeff_next*temp%next%f_cell%numerical_regularization(1)

  case('previous')
   p_cell%or_state=p_cell%or_state-0.5d0*coeff_previous*temp%previous%f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-0.5d0*coeff_current*f_cell%numerical_regularization(2)
   coeff_previous_previous=dabs(temp%previous%previous%f_cell%numer_reg_coefficient)
   p_cell%or_state=p_cell%or_state+0.5d0*coeff_previous_previous*temp%previous%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+0.5d0*coeff_next*temp%next%f_cell%numerical_regularization(1)

  case('next    ')
   p_cell%or_state=p_cell%or_state-0.5d0*coeff_current*f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-0.5d0*coeff_next*temp%next%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+0.5d0*coeff_previous*temp%previous%f_cell%numerical_regularization(2)
   coeff_next_next=dabs(temp%next%next%f_cell%numer_reg_coefficient)
   p_cell%or_state=p_cell%or_state+0.5d0*coeff_next_next*temp%next%next%f_cell%numerical_regularization(1)

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