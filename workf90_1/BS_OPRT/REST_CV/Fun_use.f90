module function_in_use
! This module provides the form of the 2D function to be integrated
! and the upper and lower limits of the integral.

use physical_state
! The file containing the module 'physical_state'.

implicit none

type(state) :: u_constant, slop_x, slop_y
real(8) :: a_upper, b_upper, c_upper, a_lower, b_lower, c_lower

public  u_constant, slop_x, slop_y, a_upper, b_upper, c_upper, &
        a_lower, b_lower, c_lower, function_in_2D, upper_limit, &
		lower_limit


contains


function function_in_2D(x,y) result(c)
! The integrand is a linear function of two variables x and y. 'u_constant' is the constant,
! 'slop_x' and 'slop_y' are the two slops in the two directions. 

implicit none
real(8), intent(in) :: x, y
type(state) :: c

c=u_constant+x*slop_x+y*slop_y

end function function_in_2D


function upper_limit(t) result(c)
! The upper limit of the integral is a second-order polynomial with 'a_upper', 'b_upper'
! and 'c_upper' as its coefficients.

implicit none
real(8), intent(in) :: t
real(8) :: c

real(8) :: cc

cc=c_upper
if(b_upper.gt.-90.0d0) then
 cc=cc+b_upper*t
 if(a_upper.gt.-90.0d0) cc=cc+a_upper*t*t
else
 if(a_upper.gt.-90.0d0) then
  print*, 'There must be something wrong!'; pause
 end if
end if

c=cc

end function upper_limit


function lower_limit(t) result(c)
! The upper limit of the integral is a second-order polynomial with 'a_lower', 'b_lower'
! and 'c_lower' as its coefficients.

implicit none
real(8), intent(in) :: t
real(8) :: c

real(8) :: cc

cc=c_lower
if(b_lower.gt.-90.0d0) then
 cc=cc+b_lower*t
 if(a_lower.gt.-90.0d0) cc=cc+a_lower*t*t
else
 if(a_lower.gt.-90.0d0) then
  print*, 'There must be something wrong!'; pause
 end if
end if

c=cc

end function lower_limit


end module function_in_use