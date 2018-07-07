module polynomials
! This module is for constructing linear and quadric polynomials.

implicit none

public  polynomial, polynomial_roots, polynomial_function

 
contains


subroutine polynomial(values,points,coefficients)
! Find interpolation linear or quadric polynomials based on given
! values at given points.

implicit none
real(8), dimension(:), intent(in) :: values, points
real(8), dimension(:), intent(out) :: coefficients

integer :: order

coefficients=-1.0d3
order=size(values)-1

select case(order)

 case(1)
  coefficients(2)=(values(2)-values(1))/(points(2)-points(1))
  coefficients(1)=(values(1)*points(2)-values(2)*points(1))/ &
                (points(2)-points(1))

 case(2)
  coefficients(3)=values(1)/(points(1)-points(2))/(points(1)- &
                  points(3))
  coefficients(3)=coefficients(3)+values(2)/(points(2)- &
                  points(1))/(points(2)-points(3))
  coefficients(3)=coefficients(3)+values(3)/(points(3)- &
                  points(1))/(points(3)-points(2))

  coefficients(2)=-values(1)*(points(2)+points(3))/(points(1)- &
                  points(2))/(points(1)-points(3))
  coefficients(2)=coefficients(2)-values(2)*(points(3)+points(1)) &
                  /(points(2)-points(1))/(points(2)-points(3))
  coefficients(2)=coefficients(2)-values(3)*(points(1)+points(2)) &
                  /(points(3)-points(1))/(points(3)-points(2))

  coefficients(1)=values(1)*points(2)*points(3)/(points(1)- &
                  points(2))/(points(1)-points(3))
  coefficients(1)=coefficients(1)+values(2)*points(3)*points(1) &
                  /(points(2)-points(1))/(points(2)-points(3))
  coefficients(1)=coefficients(1)+values(3)*points(1)*points(2) &
                  /(points(3)-points(1))/(points(3)-points(2))

end select 

end subroutine polynomial


subroutine polynomial_roots(coefficients,right_term,roots)
! Solutions to linear or quadric equations.

implicit none
real(8), dimension(:), intent(in) :: coefficients
real(8), intent(in) :: right_term
real(8), dimension(:), intent(out) :: roots
real(8), dimension(2) :: dd

real(8) :: delta
integer :: order, i

order=size(coefficients)-1

select case(order)

 case(1)
  roots(1)=(right_term-coefficients(1))/coefficients(2)

 case(2)
  delta=coefficients(2)*coefficients(2)
  delta=delta-4.d0*coefficients(3)*(coefficients(1)-right_term)
  delta=dsqrt(delta)
  dd(1)=-coefficients(2)+delta
  dd(2)=-coefficients(2)-delta
  do i=1,2
   if(dabs(dd(i)).lt.1.0d-5) dd(i)=dsign(1.0d-5,dd(i))
   roots(i)=2.0d0*(coefficients(1)-right_term)/dd(i)
  end do
end select

end subroutine polynomial_roots

function polynomial_function(coefficients,variable) result(c)

implicit none
real(8), dimension(:), intent(in) :: coefficients
real(8), intent(in) :: variable
real(8) :: c

integer :: order

order=size(coefficients)-1

c=coefficients(1)+coefficients(2)*variable
if(order.eq.2) then
 c=c+coefficients(3)*variable*variable
end if

end function polynomial_function 
  

end module polynomials