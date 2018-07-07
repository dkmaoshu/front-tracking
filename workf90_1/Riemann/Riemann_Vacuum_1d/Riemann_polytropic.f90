module Riemann_polytropic_Phi
! Phi function and computation density through a rarefaction.

use polytropic

implicit none

public  phi_polytropic, density_rare_polytropic


contains


function phi_polytropic(rhol,pl,ul,p_star) result(m)

implicit none

real(8), intent(in) :: rhol, pl
type(polytro), intent(in) :: ul
real(8), intent(in) :: p_star
real(8) :: m

real(8) :: x, gamma, epslon

x=p_star/pl
epslon=1.d-6
gamma=ul%thermo%gamma

if(x.ge.1.0d0-epslon)then
 m=dsqrt((gamma+1.d0)/2.0d0*x+(gamma-1.d0)/2.0d0)
else
 m=(gamma-1.0d0)/(2.0d0*dsqrt(gamma))*(1.0d0-x)/ &
   (1.0d0-x**((gamma-1.d0)/(2.0d0*gamma)))
end if

m=dsqrt(pl*rhol)*m

end function phi_polytropic


function density_rare_polytropic(rho_side,p_side,u_side,p_star) result(c)

implicit none

real(8), intent(in) :: p_star, rho_side, p_side
type(polytro) :: u_side
real(8) :: c

real(8) :: aa, gamma

gamma=u_side%thermo%gamma

aa=p_side/(rho_side**gamma)
c=(p_star/aa)**(1.d0/gamma)

end function density_rare_polytropic


end module Riemann_polytropic_Phi