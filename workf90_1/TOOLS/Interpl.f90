module interpolation
! This module is for interpolation of real functions.

implicit none

public  lagrange, lagrange_derivative, interpol_via_prime_fun
private lagrange_basis, lagrange_basis_derivative


contains


function lagrange(u,x,y,n) result(c)
! Lagrange interpolation function, where the array $x$ is the interpolation points, the 
! array $u$ is the values of the interpolated function at the points and $n$ is the number
! of the interpolation points. The interoplation polynomial is evaluated at $y$.

implicit none

integer, intent(in) :: n
real(8), dimension(0:n), intent(in) :: u, x
real(8), intent(in) :: y
real(8) :: c

integer :: i
real(8), dimension(0:n) :: p

c=0.0d0; p=0.0d0

do i=0,n
 p(i)=lagrange_basis(i,x,y,n)
end do
do i=0,n
 c=c+p(i)*u(i)
end do

end function lagrange


function lagrange_derivative(u,x,y,n) result(c)
! Derivative of Lagrange interpolation function, where the array $x$ 
! is the interpolation points, the array $u$ is the values of the 
! interpolated function at the points and $n$ is the number	of the 
! interpolation points. The derivative of the interoplation 
! polynomial is evaluated at $y$.

implicit none

integer, intent(in) :: n
real(8), dimension(0:n), intent(in) :: u, x
real(8), intent(in) :: y
real(8) :: c

integer :: i
real(8), dimension(0:n) :: p

c=0.0d0; p=0.0d0

do i=0,n
 p(i)=lagrange_basis_derivative(i,x,y,n)
end do
do i=0,n
 c=c+p(i)*u(i)
end do

end function lagrange_derivative 


function lagrange_basis(i,x,y,n) result(c)
! The basis functions of Lagrange interpolation.

implicit none
integer, intent(in) :: i, n
real(8), dimension(0:n), intent(in) :: x
real(8), intent(in) :: y
real(8) :: c

integer :: j

c=1.0d0
do j=0,n
 if(j.ne.i) then
  if(dabs(x(i)-x(j)).lt.0.1e-12) then
   print*, 'Something must be wrong in Interpolation'; pause
  end if
  c=c*(y-x(j))/(x(i)-x(j))
 end if
end do

end function lagrange_basis


function lagrange_basis_derivative(i,x,y,n) result(c)
! The basis functions of Lagrange interpolation.

implicit none
integer, intent(in) :: i, n
real(8), dimension(0:n), intent(in) :: x
real(8), intent(in) :: y
real(8) :: c

integer :: j, k
real(8) :: denominator, numerator
real(8), dimension(0:n) :: pp

denominator=1.0d0
do j=0,n
 if(j.ne.i) then
  if(dabs(x(i)-x(j)).lt.0.1e-12) then
   print*, 'Something must be wrong in Interpolation'; pause
  end if
  denominator=denominator*(x(i)-x(j))
 end if
end do
do k=0,n
 if(k.ne.i) then; pp(k)=1.0d0; else; pp(k)=0.0d0; end if
 do j=0,n
  if(j.ne.i.and.j.ne.k) pp(k)=pp(k)*(y-x(j))
 end do
end do
numerator=0.0d0
do k=0,n
 numerator=numerator+pp(k)
end do
c=numerator/denominator
   
end function lagrange_basis_derivative


function interpol_via_prime_fun(u,x,y,n,first_half_width,direction) result(c)
! Interpolation via prime function, where the array $x$ is the centers of interpolation cells, 
! the array $u$ is the cell-averages of the interpolated function in the cells and $n$ is the 
! number of the interpolation cells. The interoplation polynomial is evaluated at $y$.
! One more variable needed is the half length of the first grid cell, 'h'.

implicit none

integer, intent(in) :: n
real(8), dimension(0:n), intent(in) :: u, x
real(8), intent(in) :: y
real(8), intent(in) :: first_half_width
character*8, intent(in) :: direction
real(8) :: c

real(8), dimension(0:n+1) :: u_prime, x_prime
integer :: i, j
real(8) :: half_width

x_prime=0.0d0; u_prime=0.0d0
! Produce the grid points for the prime function.
half_width=first_half_width
select case(direction)

 case('upward  ')
  x_prime(0)=x(0)-half_width
  do i=1,n
   half_width=x(i)-x(i-1)-half_width
   x_prime(i)=x(i)-half_width
  end do
  x_prime(n+1)=x(n)+half_width

 case('downward')
  x_prime(n+1)=x(0)+half_width
  do i=n,1,-1
   half_width=x(n-i)-x(n+1-i)-half_width
   x_prime(i)=x(n+1-i)+half_width
  end do
  x_prime(0)=x(n)-half_width

 case default
  print*, 'Something must be wrong in interpolation!'; pause
 
end select

! Produce the data of the prime function at grid points.
select case(direction)

case('upward  ')
 u_prime=0.0d0
 do i=1,n+1; do j=1,i
  u_prime(i)=u_prime(i)+u(j-1)*(x_prime(j)-x_prime(j-1))
 end do; end do  

case('downward')
 u_prime=0.0d0
 do i=1,n+1; do j=1,i
  u_prime(i)=u_prime(i)+u(n+1-j)*(x_prime(j)-x_prime(j-1))
 end do; end do
 
end select    

c=lagrange_derivative(u_prime,x_prime,y,n+1)

end function interpol_via_prime_fun


end module interpolation