module numerical_NS_diffusion
! This module is for the Navier-Stokes numerical diffusion on the interface.

use solution
! 'solution.f90'

implicit none

public  num_NS_diffusion
!private g_cell, p_cell, l_store, r_store, cleaning_data, return_data, &
!        uniform_states_between, recleaning_data_4_stack, return_data_4_stack


contains


subroutine num_NS_diffusion
! Implement Navier-Stokes type numerical diffusion on the tracked interface.

integer :: i

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call diffusion_on_curve(i)
 end if
end do
  
do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call data_return(i)
 end if
end do

end subroutine num_NS_diffusion


subroutine diffusion_on_curve(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
integer :: head_mark
logical :: head_switch

real(8), dimension(-1:1,-1:1) :: x_velo, y_velo, pres
character*6, dimension(-1:1,-1:1) :: local_map
type(state), dimension(4) :: flv, flp
integer :: i0, j0, l
real(8) :: mu, kap

!mu=0.1d0*rr; kap=0.1d0*rr
!mu=0.25d0*rr; kap=0.25d0*rr

!mu=0.5d0*rr; kap=0.5d0*rr
mu=viscosity*rr; kap=heat_conduction*rr

nullify(temp)

if(associated(cvv(cv_nmb)%begin%next)) then
 temp=>cvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch)

 g_cell=temp%g_cell; p_cell=temp%p_cell
 i0=g_cell%x_idx; j0=g_cell%y_idx
 temp%memo=error_data

 x_velo=error_data; y_velo=error_data; pres=error_data
   
 do l=1,4
  flv(l)%gamma=error_data; flv(l)%value=0.0d0
  flp(l)%gamma=error_data; flp(l)%value=0.0d0
 end do
   
 call compute_velocity_pressure
   
! Compute diffusion fluxes.
 if(local_map(0,-1).eq.'crit  ') then  
  flv(1)%value(2)=x_velo(0,-1)-x_velo(0,0)
  flv(1)%value(3)=y_velo(0,-1)-y_velo(0,0)
  flv(1)%value(4)=0.5d0*(x_velo(0,-1)**2.0d0-x_velo(0,0)**2.0d0+y_velo(0,-1)**2.0d0-y_velo(0,0)**2.0d0)
  flp(1)%value(4)=pres(0,-1)-pres(0,0)
 end if  
  
 if(local_map(1,0).eq.'crit  ') then
  flv(2)%value(2)=x_velo(1,0)-x_velo(0,0)
  flv(2)%value(3)=y_velo(1,0)-y_velo(0,0)
  flv(2)%value(4)=0.5d0*(x_velo(1,0)**2.0d0-x_velo(0,0)**2.0d0+y_velo(1,0)**2.0d0-y_velo(0,0)**2.0d0)
  flp(2)%value(4)=pres(1,0)-pres(0,0)
 end if
   
 if(local_map(0,1).eq.'crit  ') then
  flv(3)%value(2)=x_velo(0,1)-x_velo(0,0)
  flv(3)%value(3)=y_velo(0,1)-y_velo(0,0)
  flv(3)%value(4)=0.5d0*(x_velo(0,1)**2.0d0-x_velo(0,0)**2.0d0+y_velo(0,1)**2.0d0-y_velo(0,0)**2.0d0)
  flp(3)%value(4)=pres(0,1)-pres(0,0)
 end if
  
 if(local_map(-1,0).eq.'crit  ') then    
  flv(4)%value(2)=x_velo(-1,0)-x_velo(0,0)
  flv(4)%value(3)=y_velo(-1,0)-y_velo(0,0)
  flv(4)%value(4)=0.5d0*(x_velo(-1,0)**2.0d0-x_velo(0,0)**2.0d0+y_velo(-1,0)**2.0d0-y_velo(0,0)**2.0d0)
  flp(4)%value(4)=pres(-1,0)-pres(0,0)
 end if
  
 temp%memo=p_cell%or_state
     
 do l=1, 4
  temp%memo=temp%memo+mu*flv(l)+kap*flp(l)
 end do
        
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do


contains


subroutine compute_velocity_pressure

implicit none

type(critical_cell), pointer :: temp_g

integer:: ii, jj

x_velo=error_data; y_velo=error_data; pres=error_data
local_map='bound '

do ii=-1, 1
 do jj=-1, 1
  if(i0+ii.ge.nxll.and.i0+ii.le.nxll+nx-1) then
   if(j0+jj.ge.nyll.and.j0+jj.le.nyll+ny-1) then
      
	select case(ggd_cell(i0+ii,j0+jj)%region) 
      
     case('smth  ')
      local_map(ii,jj)='smth  '
      
     case('crit  ') 
	  local_map(ii,jj)='crit  '
      call visit(ggd_cell(i0+ii,j0+jj)%ccpt(1)%address,temp_g)	  	     
      x_velo(ii,jj)=x_velocity(temp_g%p_cell%l_state)
      y_velo(ii,jj)=y_velocity(temp_g%p_cell%l_state)
      pres(ii,jj)=pressure(temp_g%p_cell%l_state)    
        
     case default; call error_message
      	  
    end select 	   
     
   end if
  end if
 end do
end do 

end subroutine compute_velocity_pressure


end subroutine diffusion_on_curve


subroutine data_return(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
integer :: head_mark
logical :: head_switch

nullify(temp)

if(associated(cvv(cv_nmb)%begin%next)) then
 temp=>cvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch)

 g_cell=temp%g_cell; p_cell=temp%p_cell
  
 temp%p_cell%or_state=temp%memo
 temp%memo=error_data
     
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine data_return


end module numerical_NS_diffusion