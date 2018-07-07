module tools
! This module provides some basic mathematical tools, such as
! finding the intersection point of two lines, compute area of 
! polygons and so on.

use polynomials
! 'polynml.f90'

implicit none

public  ints_of_lines, area_of_polygon, length_of_segment, &
        normal_of_line, line_pass_points, cycle_4, pivots_deletion, between
private area_of_triangle


contains


subroutine ints_of_lines(pts1,pts2,pt,error_message)
! This subroutine determins the intersection point of two lines 
! given by two points on each of them. In this function, pts1(1,:) 
! and pts1(2,:), and pts2(1,:) and pts2(2,:) are the two pairs of 
! points on the first and second lines and pt is the intersection 
! point.

implicit none
real(8), dimension(2,2), intent(in) :: pts1, pts2
real(8), dimension(2), intent(out) :: pt
character*10, intent(out) :: error_message

real(8) :: a11, a12, a21, a22, c1, c2, dd, r1, r2

error_message='ok        '; pt=-1000.0d0
r1=(pts1(1,1)-pts1(2,1))*(pts1(1,1)-pts1(2,1))+(pts1(1,2)-pts1(2,2))*(pts1(1,2)-pts1(2,2))
r1=dsqrt(r1)
r2=(pts2(1,1)-pts2(2,1))*(pts2(1,1)-pts2(2,1))+(pts2(1,2)-pts2(2,2))*(pts2(1,2)-pts2(2,2))
r2=dsqrt(r2)
if(r1.lt.1.0d-2.or.r2.lt.1.0d-2) then
 error_message='short     '; return
end if

a11=-(pts1(1,2)-pts1(2,2))
a12=(pts1(1,1)-pts1(2,1))
c1=(pts1(1,1)-pts1(2,1))*pts1(1,2)-(pts1(1,2)-pts1(2,2))*pts1(1,1)
a21=-(pts2(1,2)-pts2(2,2))
a22=(pts2(1,1)-pts2(2,1))
c2=(pts2(1,1)-pts2(2,1))*pts2(1,2)-(pts2(1,2)-pts2(2,2))*pts2(1,1)
dd=a11*a22-a12*a21
if(dabs(dd).lt.0.1d-4) then
 error_message='parallel  '
 return
end if 
pt(1)=(c1*a22-c2*a12)/dd
pt(2)=(a11*c2-a21*c1)/dd
return      

end	subroutine ints_of_lines


function area_of_triangle(pt1,pt2,pt3)
! This function calculates the area of a triangular given by 
! its three vortices. In this function, pt1, pt2 and pt3 are 
! the three vortices.      

implicit none
real(8), dimension(2), intent(in) :: pt1, pt2, pt3
real(8) :: area_of_triangle

area_of_triangle=(pt1(1)*pt2(2)-pt2(1)*pt1(2))+ &
                 (pt2(1)*pt3(2)-pt3(1)*pt2(2))+ &
				 (pt3(1)*pt1(2)-pt1(1)*pt3(2))
area_of_triangle=0.5d0*dabs(area_of_triangle)                

end	function area_of_triangle


function area_of_polygon(pt,n) result(c)
! This function calculates the area of a polygon given by 
! its 'n' vortices. In this function, 'pt' are the vortices and 
! 'n' is the number of the vertices of the polygon.      
 
implicit none
integer, intent(in) :: n
real(8), dimension(n,2), intent(in) :: pt
real(8) :: c

real(8), dimension(2) :: pt1, pt2, pt3 
integer :: i

c=0.d0
do i=2,n-1
 pt1=pt(1,:); pt2=pt(i,:); pt3=pt(i+1,:)
 c=c+area_of_triangle(pt1,pt2,pt3)
end do

end function area_of_polygon       


function length_of_segment(pt1,pt2) result(c)
! This function computes the length of a line segment given by 
! its two end points 'pt1' and 'pt2'.

real(8), dimension(2), intent(in) :: pt1, pt2
real(8) :: c

c=(pt2(1)-pt1(1))*(pt2(1)-pt1(1))+(pt2(2)-pt1(2))*(pt2(2)-pt1(2))
c=dsqrt(c)

end function length_of_segment


function normal_of_line(pt1,pt2) result(c)
! This function computes the normal of a line given by two points
! on it. In this function, 'pt1' is the first point and 'pt2' is 
! the second point.

implicit none

real(8), dimension(2), intent(in) :: pt1, pt2
real(8), dimension(2) :: c

real(8) :: length

length=length_of_segment(pt1,pt2)

c(1)=(pt2(2)-pt1(2))/length
c(2)=-(pt2(1)-pt1(1))/length

end function normal_of_line


function cycle_4(k)

implicit none
integer :: k, cycle_4

if(k.lt.1) then
 cycle_4=k+4
else
 if(k.gt.4) then
  cycle_4=k-4
 else
  cycle_4=k
 end if
end if

end function cycle_4


function distance_to_line(a,b,c,pt) 
! This function computs the distance of a given point to a line
! of 'ax+by+c=0'.

implicit none
real(8), intent(in)	:: a, b, c
! The coefficients of the line equation.
real(8), dimension(2) :: pt
! The input point.
real(8) :: distance_to_line

distance_to_line=a*pt(1)+b*pt(2)+c

end function distance_to_line


subroutine line_pass_points(pt1,pt2,a,b,c)
! The subroutine computes the coefficients of the line pass two
! given points.

implicit none
real(8), dimension(2), intent(in) :: pt1, pt2
! The two passed points
real(8), intent(out) :: a, b, c
! The coefficients of the line equation, which takes the form
! 'ax+by+c=0'.

a=pt1(2)-pt2(2)
b=pt2(1)-pt1(1)
c=pt1(1)*pt2(2)-pt2(1)*pt1(2)

end subroutine line_pass_points


subroutine pivots_deletion(q,bb,x,sma)

implicit none

real(8), dimension(:,:), intent(in) ::  q
real(8), dimension(:), intent(in) :: bb
real(8), dimension(:), intent(out) :: x, sma

real(8), dimension(:,:), allocatable :: qq
real(8), dimension(:), allocatable :: b
real(8), dimension(:), allocatable :: y
integer, dimension(:), allocatable :: llr
integer :: k, l, ls, kr, lr, ll, size_of_unknown
real(8) :: xx

size_of_unknown=size(bb)
allocate(qq(size_of_unknown,size_of_unknown))
allocate(b(size_of_unknown))
allocate(y(size_of_unknown))
allocate(llr(size_of_unknown))

do l=1,size_of_unknown
 x(l)=0.0d0
 y(l)=0.0d0
 sma(l)=0.0d0
 llr(l)=l
 b(l)=bb(l)
end do

do l=1,size_of_unknown
 do k=1,size_of_unknown
  qq(l,k)=q(l,k)
 end do
end do

!do nk1=1,3
! write(*,'(5f14.4)') qq(nk1,1),qq(nk1,2),qq(nk1,3),qq(nk1,4),b(nk1)
!end do
!print*

do ls=1,size_of_unknown
 do l=ls,size_of_unknown
  do k=ls,size_of_unknown
   xx=dabs(qq(k,l))
   if(sma(ls).lt.xx) then
    sma(ls)=xx
    kr=k
    lr=l
   endif
  end do
 end do
 if(sma(ls).gt.sma(1)*0.1d-6) then
  do l=ls,size_of_unknown
   xx=qq(kr,l)
   qq(kr,l)=qq(ls,l)
   qq(ls,l)=xx
  end do
  xx=b(kr)
  b(kr)=b(ls)
  b(ls)=xx
  do k=1,size_of_unknown
   xx=qq(k,lr)
   qq(k,lr)=qq(k,ls)
   qq(k,ls)=xx
  end do
  ll=llr(ls)
  llr(ls)=llr(lr)
  llr(lr)=ll
  do k=ls+1,size_of_unknown
   xx=qq(k,ls)/qq(ls,ls)
   do l=ls,size_of_unknown
    qq(k,l)=qq(k,l)-xx*qq(ls,l)
   end do
   b(k)=b(k)-xx*b(ls)
  end do
 else
  exit
 endif

! do nk1=1,3
!  write(*,'(5f14.4)') qq(nk1,1),qq(nk1,2),qq(nk1,3),qq(nk1,4),b(nk1)
! end do
! print*

end do

do k=ls-1,1,-1
 do ll=k+1,ls-1
  y(k)=y(k)-y(ll)*qq(k,ll)
 end do
 y(k)=(y(k)+b(k))/qq(k,k)
end do
do l=1,ls-1
 do k=1,size_of_unknown
  x(llr(l))=y(l)
 end do
end do
do l=1,size_of_unknown
 sma(l)=1.0d0
end do
do l=ls,size_of_unknown
 sma(llr(l))=0.
end do

deallocate(qq,b,y,llr)

end subroutine pivots_deletion


function between(x,a,b) result(c)

implicit none

real(8), intent(in) :: x, a, b
real(8) :: c

real(8) :: aa, bb

aa=dmin1(a,b); bb=dmax1(a,b)

c=x; if(x.lt.aa) c=aa; if(x.gt.bb) c=bb

end function between


end module tools