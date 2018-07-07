module numerical_fluxes
! This module is for the numerical fluxes for either scalar 
! conservation law or Euler system. The type-check of parameter is 
! done by interface

use physical_state
use flux_functions
implicit none

type(state) :: maximum_characteristic_speeds

public flux
private wflux_x, wflux_y


contains


function flux(u,direction,x_posi,y_posi) result(c)

! Function flux computes the numerical flux on given physical_states
! 'u(1)', 'u(2)', 'u(3)' and 'u(4)'. The argument 'direction' 
! decides whether the flux in x- or y-direction is computed. WENO flux 
! is used for the numerical flux, rf. "Efficient implementation of weighted
! ENO schemes", by Guang-Shan Jiang and Chi-Wang Shu, in Journal of 
! Computational Physics, v126 (1996), pp.202-228.

implicit none
type(state), intent(in), dimension(-2:3) :: u
character*2, intent(in) :: direction
real(8), intent(in) :: x_posi, y_posi
type(state) :: c

real(8), dimension(-2:3,4) :: wu
real(8), dimension(4) :: wf
integer :: i

emxy=maximum_characteristic_speeds%value

do i=-2, 3
 wu(i,:)=u(i)%value
end do
select case(direction)
 case('xx')
  call wflux_x(wu,wf,u(0)%gamma)
 case('yy')
  call wflux_y(wu,wf,u(0)%gamma)
 case default; call error_message
end select

c%value=wf
c%gamma=u(0)%gamma

end function flux

end module numerical_fluxes