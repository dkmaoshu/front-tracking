module flux_functions
! This module provides all the functions used in the construction of
! the second-odrer TVD numerical flux function.

use grid
use euler_functions

implicit none

public tflux


contains


subroutine tflux(w4,w3,w2,w1,flux)

! The subroutine tflux constructs a TVB flux by a Roe flux. The way 
! the flux is constructed is suggested by Osher (cf.Osher's "Very 
! high order accurate TVD scheme"). The TVB improvement is suggested
! by Shu.

implicit none

real(8), dimension(4), intent(in) :: w1, w2, w3, w4
real(8), dimension(4), intent(out) :: flux 

real(8), dimension(4) :: char1, coef1, char2, coef2, dfn1, dfn0,  &
                         dfp0, dfpn1, aa3, aa4, he
real(8), dimension(4,4) :: e1, e2, el1, el2
real(8), parameter :: alpha1=0.1666666d0, b=2.0d0

integer :: k, l
real(8) :: aa1, aa2

dfn1=0.d0
dfn0=0.d0
dfp0=0.d0
dfpn1=0.d0
call roe(w4,w3,el1,e1,char1,coef1)
call roe(w3,w2,el2,e2,char2,coef2)

do k=1,4
 aa1=snega(char1(k))*coef1(k)
 aa2=snega(char2(k))*coef2(k)
 aa3(k)=amid(aa1,b*aa2)
 aa4(k)=amid(aa2,b*aa1)
end do  
do l=1,4
 do k=1,4
  dfn1(l)=dfn1(l)+aa3(k)*e1(k,l)
  dfn0(l)=dfn0(l)+aa4(k)*e2(k,l)
 end do  
end do  

call roe(w2,w1,el1,e1,char1,coef1)
do k=1,4
 aa2=sposi(char2(k))*coef2(k)
 aa1=sposi(char1(k))*coef1(k)
 aa3(k)=amid(aa2,b*aa1)
 aa4(k)=amid(aa1,b*aa2)
end do
do l=1,4
 do k=1,4
  dfp0(l)=dfp0(l)+aa3(k)*e2(k,l)
  dfpn1(l)=dfpn1(l)+aa4(k)*e1(k,l)
 end do  
end do

call rflux(w3,w2)
flux=he-alpha1*dfn1-(0.5d0-alpha1)*dfn0+(0.5d0-alpha1)*dfp0+ &
     alpha1*dfpn1


contains


subroutine rflux(w2,w1)
! The subroutine rflux constructs a first Roe scheme flux.

implicit none

real(8), dimension(4), intent(in) :: w1, w2

real(8), dimension(4,4) :: e, el
real(8), dimension(4) :: ch, coef
integer :: k

call roe(w2,w1,el,e,ch,coef)
do k=1,4
 he(k)=0.5d0*(flux_in_x(k,w2)+flux_in_x(k,w1))
 do l=1,4
  he(k)=he(k)-0.5d0*coef(l)*dabs(ch(l))*e(l,k)
 end do
end do

end subroutine rflux


function amid(x,y)
implicit none
real(8), intent(in) :: x, y
real(8) :: amid
real(8) :: z, yy
!yy=y+50.d0*dsign(16.0d0*h*h,x)
yy=y
if(dabs(x).le.dabs(yy)) then
 z=x
else
 z=yy
end if
if(x*yy.gt.0.0d0) then
 amid=z
else
 amid=0.d0
end if
end function amid


function snega(x)
implicit none
real(8), intent(in) :: x
real(8) :: snega
if(x.le.0.0d0) then
 snega=x
else
 snega=0.0d0
end if
end	function snega


function sposi(x)
implicit none
real(8), intent(in) :: x
real(8) :: sposi
if(x.le.0.0d0) then
 sposi=0.0d0
else
 sposi=x
end if
end	function sposi


end subroutine tflux


end module flux_functions