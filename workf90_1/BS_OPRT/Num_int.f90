module numerical_integrals
! This module provides several formulas of numerical integrals

use physical_state
! 'linear.f90'

implicit none


contains


subroutine numerical_integration(f,a,b,method,value)

implicit none
type(state), external :: f
! The function to be integrated.
real(8), intent(in) :: a ,b
! The interval of integral.
character*10, intent(in) :: method
! The numerical method  used.
type(state), intent(out) :: value
! The value of the numerical integral.

type(state) :: fff
real(8) :: ab

fff=-1000.0d0; value=0.0d0

select case(method)

 case('trapezoid ')
  fff=f(a)
  value=value+0.5d0*(b-a)*fff
  fff=f(b)
  value=value+0.5d0*(b-a)*fff

 case('simpson   ')
  fff=f(a)
  value=value+(b-a)/6.0d0*fff
  fff=f(b)
  value=value+(b-a)/6.0d0*fff
  ab=0.5d0*(a+b)
  fff=f(ab)
  value=value+4.0d0*(b-a)/6.0d0*fff
 
 case default
  print*, 'There must be something wrong!'; pause

end select

end subroutine numerical_integration


end module numerical_integrals