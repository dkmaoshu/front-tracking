module numerical_fluxes
! This module is for the numerical fluxes for either scalar 
! conservation law or Euler system. The type-check of parameter is 
! done by interface

use physical_state
use flux_functions
implicit none

public flux
private tflux


contains


function flux(u,direction,x_posi,y_posi) result(c)

! Function flux computes the numerical flux on given physical_states
! 'u(1)', 'u(2)', 'u(3)' and 'u(4)'. The argument 'direction' 
! decides whether the flux in x- or y-direction is computed. The 
! functions he and amid are related functions used in the flux 
! functions. ! Refer to Osher & Chakravarthy's `High resolution 
! schemes and entropy condition', SIAM J. Numer. Anal. 21 (1984), 
! pp. 955-984 for detail description of the method.

implicit none
type(state), intent(in), dimension(4) :: u
character*2, intent(in) :: direction
real(8), intent(in) :: x_posi, y_posi
type(state) :: c

real(8), dimension(4) :: w1, w2, w3, w4, ww

w1=u(1)%value; w2=u(2)%value; w3=u(3)%value; w4=u(4)%value
select case(direction)
 case('xx')
  call tflux(w4,w3,w2,w1,ww)
 case('yy')
  call shift(w1); call shift(w2); call shift(w3); call shift(w4)
  call tflux(w4,w3,w2,w1,ww)
  call shift(ww)
 case default; call error_message
end select

c%value=ww


contains


subroutine shift(w)
implicit none
real(8), dimension(4), intent(inout) :: w
real(8) :: xx

xx=w(2); w(2)=w(3); w(3)=xx

end subroutine shift


end function flux

end module numerical_fluxes