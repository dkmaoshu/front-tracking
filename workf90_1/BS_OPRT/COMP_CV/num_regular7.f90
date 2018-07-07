module numerical_regularization
! This module is for the numerical regularization on discontinuity curves.

! Version 7: No mass difussion. Diffusions occur only in the momentum and total energy to dissipate the velocity and pressure.

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

real(8) :: l_area, r_area, l_area_p, r_area_p, length, length_p
type(state) :: physical_difference, or_state_temp, or_state_temp_p, ttt, tttl, tttr

! Code for version 7. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
real(8) :: uu, vv, pp, uu_p, vv_p, pp_p, ddensity, ddensity_p, ccoef, ccoef_p
real(8) :: dd1, dd2, ddp1, ddp2, alpha_l, alpha_r
    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
 p_cell=temp%p_cell; g_cell=temp%g_cell
  
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
    
! Code for version 7. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
   uu=l_area*x_velocity(p_cell%l_state)+r_area*x_velocity(p_cell%r_state)
   vv=l_area*y_velocity(p_cell%l_state)+r_area*y_velocity(p_cell%r_state)
   pp=l_area*pressure(p_cell%l_state)+r_area*pressure(p_cell%r_state)
   
   call obtain_or(p_cell%or_state%value(1),uu,vv,pp,p_cell%l_state%gamma,p_cell%r_state%gamma,l_area,r_area,or_state_temp)
    
!   call coefficients(p_cell%l_state%value(1),p_cell%r_state%value(1),p_cell%or_state%value(1),alpha_l,alpha_r) 
    
!   call tran3(p_cell%l_state%value(1),uu,vv,pp,tttl%value,p_cell%l_state%gamma)
!   call tran3(p_cell%r_state%value(1),uu,vv,pp,tttr%value,p_cell%r_state%gamma)
!   or_state_temp=alpha_l*tttl+alpha_r*tttr
    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       
!   or_state_temp=l_area*p_cell%l_state+r_area*p_cell%r_state
   physical_difference=p_cell%or_state-or_state_temp
   
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
    
! Code for version 7. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   uu=l_area*x_velocity(p_cell%l_state)+r_area*x_velocity(p_cell%r_state)
   vv=l_area*y_velocity(p_cell%l_state)+r_area*y_velocity(p_cell%r_state)
   pp=l_area*pressure(p_cell%l_state)+r_area*pressure(p_cell%r_state)
   
   uu_p=l_area_p*x_velocity(p_cell_p%l_state)+r_area_p*x_velocity(p_cell_p%r_state)
   vv_p=l_area_p*y_velocity(p_cell_p%l_state)+r_area_p*y_velocity(p_cell_p%r_state)
   pp_p=l_area_p*pressure(p_cell_p%l_state)+r_area_p*pressure(p_cell_p%r_state)
    
   ccoef=length_p/(length+length_p); ccoef_p=length/(length+length_p)

   dd1=l_area*p_cell%l_state%value(1)+r_area*p_cell%r_state%value(1) 
   ddp1=l_area_p*p_cell_p%l_state%value(1)+r_area_p*p_cell_p%r_state%value(1)  
   dd2=ttt%value(1)-ddp1
   ddp2=ttt%value(1)-dd1
    
   ddensity=ccoef*dd1+ccoef_p*dd2; ddensity_p=ccoef*ddp2+ccoef_p*ddp1
   
   call obtain_or(ddensity,uu,vv,pp,p_cell%l_state%gamma,p_cell%r_state%gamma,l_area,r_area,or_state_temp)
   call obtain_or(ddensity_p,uu_p,vv_p,pp_p,p_cell_p%l_state%gamma,p_cell_p%r_state%gamma,l_area_p,r_area_p,or_state_temp_p)
     
!   call coefficients(p_cell%l_state%value(1),p_cell%r_state%value(1),ddensity,alpha_l,alpha_r) 
!   call tran3(p_cell%l_state%value(1),uu,vv,pp,tttl%value,p_cell%l_state%gamma)
!   call tran3(p_cell%r_state%value(1),uu,vv,pp,tttr%value,p_cell%r_state%gamma)
!   or_state_temp=alpha_l*tttl+alpha_r*tttr
    
!   call coefficients(p_cell_p%l_state%value(1),p_cell_p%r_state%value(1),ddensity_p,alpha_l,alpha_r) 
!   call tran3(p_cell_p%l_state%value(1),uu_p,vv_P,pp_p,tttl%value,p_cell_p%l_state%gamma)
!   call tran3(p_cell_p%r_state%value(1),uu_p,vv_p,pp_p,tttr%value,p_cell_p%r_state%gamma)
!   or_state_temp_p=alpha_l*tttl+alpha_r*tttr

   physical_difference=ttt-or_state_temp-or_state_temp_p
   
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


contains


subroutine coefficients(density_l,density_r,density,coeff_l,coeff_r)

implicit none

real(8), intent(in) :: density_l, density_r, density
real(8), intent(out) :: coeff_l, coeff_r

coeff_l=(density-density_r)/(density_l-density_r)
if(coeff_l.lt.0.0d0) coeff_l=0.0d0
if(coeff_l.gt.1.0d0) coeff_l=1.0d0
coeff_r=1.0d0-coeff_l

end subroutine coefficients


subroutine obtain_or(rho,uu,vv,pp,gam1,gam2,area_1,area_2,or_statee)

implicit none

real(8), intent(in) :: rho, uu, vv, pp, gam1, gam2, area_1, area_2
type(state), intent(out) :: or_statee

or_statee%gamma=error_data

or_statee%value(1)=rho
or_statee%value(2)=rho*uu
or_statee%value(3)=rho*vv
or_statee%value(4)=0.5d0*rho*(uu*uu+vv*vv)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Only for gamma-law gas, modifications are needed for more general gases.
or_statee%value(4)=or_statee%value(4)+pp/(gam1-1.0d0)*area_1+pp/(gam2-1.0d0)*area_2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end subroutine obtain_or



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
real(8) :: coeff

!coeff=0.0d0
!coeff=rr
coeff=0.5d0*rr

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

 select case(temp%partner)

  case('single  ')
   p_cell%or_state=p_cell%or_state-coeff*f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-coeff*f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+coeff*temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+coeff*temp%next%f_cell%numerical_regularization(1)

  case('previous')
   p_cell%or_state=p_cell%or_state-coeff*temp%previous%f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-coeff*f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+coeff*temp%previous%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+coeff*temp%next%f_cell%numerical_regularization(1)

  case('next    ')
   p_cell%or_state=p_cell%or_state-coeff*f_cell%numerical_regularization(1)
   p_cell%or_state=p_cell%or_state-coeff*temp%next%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+coeff*temp%previous%f_cell%numerical_regularization(2)
   p_cell%or_state=p_cell%or_state+coeff*temp%next%next%f_cell%numerical_regularization(1)

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