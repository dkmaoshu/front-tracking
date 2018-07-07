module original_flux

use grid
implicit none
public flux_or

interface flux_or
 module procedure flux_or
end interface


contains


function flux_or(x,fg) 

implicit none

real (8) :: flux_or

real(8), parameter :: alpha1=0.1666666d0, b=2.0d0
real(8), intent(in) :: x(4) 
real(8) :: dfmp1,dfmp3,dfpm1,dfpp1,aa,a1,a2,a3,a4
real(8) :: fg

dfmp3=he(x(4),x(3),fg)-fg(x(3))
aa=he(x(3),x(2),fg)
dfmp1=aa-fg(x(2))
dfpp1=fg(x(3))-he(x(3),x(2),fg)
dfpm1=fg(x(2))-he(x(2),x(1),fg)
a4=amid(dfmp3,b*dfmp1)
a3=amid(dfmp1,b*dfmp3)
a2=amid(dfpp1,b*dfpm1)
a1=amid(dfpm1,b*dfpp1)
flux_or=he(x(3),x(2),fg)- &
alpha1*a4-(0.5-alpha1)*a3+(0.5-alpha1)*a2+alpha1*a1
	
end	function flux_or


function he(x,y,fg)

implicit none
 
real(8), intent(in) :: x,y
real(8) :: he, fg

he=0.5*(fg(x)+fg(y))
he=he-0.6*(x-y)

end	function he


function amid(x,y)

implicit none
	
real(8), intent(in) :: x,y
real(8) :: z,yy,amid

yy=y+50.*dsign(h*h,x)

if(dabs(x).le.dabs(yy)) then
  z=x
else
  z=yy
end if
if(x*yy.gt.0.) then
  amid=z
else
  amid=0.d0
end if

end	function amid

end module original_flux


! This module is for the numerical fluxes for either scalar 
! conservation law or Euler system. The type-check of parameter is 
! done by interface

module numerical_fluxes

use grid
use physical_state
use original_flux
implicit none

public flux
private flux_or, amid, he


contains


function flux(u,direction)

! Function flux computes the numerical flux on given physical_states
! 'u(1)', 'u(2)', 'u(3)' and 'u(4)'. The argument 'direction' 
! decides whether the flux in x- or y-direction is computed. The 
! functions he and amid are related functions used in the flux 
! functions. ! Refer to Osher & Chakravarthy's `High resolution 
! schemes and entropy condition', SIAM J. Numer. Anal. 21 (1984), 
! pp. 955-984 for detail description of the method.

implicit none
type(state) :: flux
type(state), dimension(4) :: u
character*2 :: direction
real(8), dimension(4) :: x(4)

x=u%value
if(direction.eq.'xx') then
 flux%value=flux_or(x,f)
else
 if(direction.eq.'yy') then
  flux%value=flux_or(x,g)
 else
  print*, 'There must be something wrong in flux.'
  stop
 end if
end if

end function flux

end module numerical_fluxes