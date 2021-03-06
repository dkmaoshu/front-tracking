module polytropic
! Polytropic gas.

use conserved_quantities
use my_messages

implicit none

type pol_thermo
 real(8) :: gamma
end type pol_thermo

type polytro
 type(conserved) :: csq
 type(pol_thermo) :: thermo
end type polytro

interface df
 module procedure density
end interface

interface mf
 module procedure momentum
end interface

interface ef
 module procedure energy
end interface

interface ief
 module procedure internal_energy
end interface

interface pf
 module procedure pressure
end interface

interface cf
 module procedure sound_speed
end interface
  
interface vf
 module procedure velocity
end interface

interface hf
 module procedure enthalpy
end interface

interface sf
 module procedure entropy
end interface


public  df, mf, ef, pf, cf, vf, hf, sf, ief
private density, momentum, energy, pressure, velocity, sound_speed, &
        entropy, enthalpy, internal_energy


contains


! The density taken as the first component of the conserved quantities.
function density(u) result(c)

implicit none

type(polytro), intent(in) :: u
real(8) :: c

c=u%csq%value(1)

end function density


! The momentum taken as the second component of the conserved quantities.
function momentum(u) result(c)

implicit none

type(polytro), intent(in) :: u
real(8) :: c

c=u%csq%value(2)

end function momentum


! The total energy taken as the third component of the conserved quantities.
function energy(u) result(c)

implicit none

type(polytro), intent(in) :: u
real(8) :: c

c=u%csq%value(3)

end function energy


! The pressure calculated from the conservative quantities. 
function pressure(u) result(c)

implicit none

type(polytro), intent(in) :: u
real(8) :: c

real(8) :: gamma

gamma=u%thermo%gamma

c=(gamma-1.0d0)*(u%csq%value(3)-0.5d0*u%csq%value(2)**2.d0/u%csq%value(1))

end function pressure


! The sound speed computed from the conservative quantities.
function sound_speed(u) result(c)

implicit none

type(polytro), intent(in):: u
real(8) :: c

real*8:: gamma

gamma=u%thermo%gamma

c=dsqrt(gamma*pf(u)/u%csq%value(1)) 

end function sound_speed


! The enthalpy calculated from the conservative quantities.
function enthalpy(u) result(c)

implicit none

type(polytro), intent(in):: u
real*8:: c
   
c=(u%csq%value(3)+pf(u))/u%csq%value(1)

end function enthalpy


! The velocity calculated from the conservative quantities.
function velocity(u) result(c)

implicit none

type(polytro), intent(in):: u
real*8:: c
 
c=u%csq%value(2)/u%csq%value(1)

end function velocity
 

! The entropy computed from the conservative quantities.
function entropy(u) result(c)
implicit none
type(polytro), intent(in):: u
real(8) :: c

real*8:: gamma


gamma=u%thermo%gamma 

c=dlog(pf(u)/u%csq%value(1)**gamma)  

end function entropy


! The internal energy computed from the conservative quantities.
function internal_energy(u) result(c)
implicit none
type(polytro), intent(in):: u
real*8:: c

c=u%csq%value(3)-0.5d0*u%csq%value(2)**2.0d0/u%csq%value(1)

end function internal_energy


! Tran3 calculates a conserved array from given density, velocity, 
! pressure, gamma, b and rho0.
subroutine dynamic_2_conserved_polytropic(rho,v,p,u)  

implicit none

real*8, intent(in) :: rho, v, p
type(polytro), intent(out) :: u 

real(8) :: gamma

gamma=u%thermo%gamma

u%csq%value(1)=rho
u%csq%value(2)=v*rho
u%csq%value(3)=p/(gamma-1.0d0)+0.5d0*rho*v**2.0d0

end subroutine dynamic_2_conserved_polytropic


end module polytropic