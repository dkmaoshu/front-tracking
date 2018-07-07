module numerical_integral_in_2D
! This module is for the 2D numerical integrations used in split of 
! critical cells. The functions used in this module are provided by
! the module 'function_in_use'.

use physical_state
! 'linear.f90'

use numerical_integrals
! 'num_int.f90'

use function_in_use
! 'fun_use.f90'

implicit none

character*3 :: geometrical_type
real(8) :: upper, lower
real(8) :: xy


contains


subroutine integral_2(value)
! The second integral.

implicit none
type(state), intent(out) :: value

call numerical_integration(integral_1,lower,upper,'simpson   ',value)

end subroutine integral_2


function integral_1(tt) result (c)
! The first integral, which produces a function with respect to
! the second variable.

implicit none
real(8), intent(in) :: tt
type(state) :: c

real(8) :: a, b
type(state) :: cc

xy=tt
a=lower_limit(tt)
b=upper_limit(tt)

call numerical_integration(reduced_function_in_1D,a,b,'trapezoid ',cc)
c=cc

end function integral_1


function reduced_function_in_1D(t) result(c)

implicit none
real(8), intent(in) :: t
type(state) :: c

select case(geometrical_type)
 case('xx ')
  c=function_in_2D(t,xy)
 case('yy ')
  c=function_in_2D(xy,t)
 case default
  print*, 'There must be something wrong!'; pause
end select

end function reduced_function_in_1D


end module numerical_integral_in_2D