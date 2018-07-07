Module 	data_extrapolation
! This module contains the tools for data extrapolation.

use grid
! 'grid.f90'

implicit none

interface data_fill
 module procedure data_fill_r
end interface

public data_fill
private data_extrap, extrapolation, data_fill_r


contains


subroutine data_fill_r(uxy,highest_order)

implicit none

real(8), dimension(-3:3), intent(inout) :: uxy
! The prepared data section.
integer, intent(in) :: highest_order
! The highest order of extrapolation.

call data_extrap(uxy,highest_order)

end subroutine data_fill_r


subroutine data_extrap(uxy,highest_order)
! This subroutine supplments the missing data in the data section
! 'uxy' with extrapolation data.

implicit none

real(8), dimension(-3:3), intent(inout) :: uxy
! The prepared data section.
integer, intent(in) :: highest_order
! The highest order of extrapolation.

integer :: l, l1, l2, l3, start_number

do l=0,3
 if(uxy(-l).gt.0.9d0*error_data) then
  start_number=-l; exit
 end if
 if(uxy(l).gt.0.9d0*error_data) then
  start_number=l; exit
 end if
end do

do l=start_number-1,-3,-1
 if(uxy(l).lt.0.9d0*error_data) exit
end do
l1=l+1
do l=start_number,3
 if(uxy(l).lt.0.9d0*error_data) exit
end do
l2=l-1
l3=l2-l1+1
if(l3.ge.highest_order) l3=highest_order
do l=-3,l1-1
 uxy(l)=extrapolation(uxy(l1),uxy(l1+1),uxy(l1+2),l1-l,l3)
end do
do l=l2+1,3
 uxy(l)=extrapolation(uxy(l2),uxy(l2-1),uxy(l2-2),l-l2,l3)
end do


end subroutine data_extrap


function extrapolation(u1,u2,u3,md,nd)
!This function computes zero, first and second extrapolation data 
! from given data 'u1', 'u2', 'u3', where 'md' the (grid) location
! of the extrapolation data and 'nd' is the order of the 
! extrapolation.

implicit none
real(8) :: vv(3,3); real(8) :: extrapolation, u1, u2, u3
integer :: md, nd

integer :: i
real(8) :: d, dd, ratio
real(8), dimension(3) :: ww

vv(1,3)=3.d0*u1-3.d0*u2+u3
vv(2,3)=6.d0*u1-8.d0*u2+3.d0*u3
vv(3,3)=10.d0*u1-15.d0*u2+6.d0*u3
!vv(4,3)=15.0d0*u1-24.0d0*u2+10.d0*u3
! Second order extrapolation.

vv(1,2)=2.d0*u1-u2
vv(2,2)=3.d0*u1-2.d0*u2
vv(3,2)=4.d0*u1-3.*u2
!vv(4,2)=5.0d0*u1-4.0d0*u2
! First order extrapolation.

vv(1,1)=u1
vv(2,1)=u1
vv(3,1)=u1
!vv(4,1)=u1
! Zero order extrapolation.

! The following is to modify the third-order extrapolation to prevent the possible 
! overshooting or collapse of the extrapolation data.
d=dabs(u2-u3); dd=dabs(u1-2.0d0*u2+u3)
do i=1, 3
 ww(i)=(2.0d0+dfloat(i))*u2-(1.0d0+dfloat(i))*u3
end do
if(d.gt.1.0d-12) then
 ratio=0.5d0*dd/(d+dd)
else
 ratio=0.5d0
end if
do i=1, 3
 vv(i,3)=(1.0d0-ratio)*vv(i,3)+ratio*ww(i)
end do

if(md.gt.3) then
 print*, 'Something is wrong in Data-Extrapolation'; pause
end if

extrapolation=vv(md,nd)

end	function extrapolation


end module data_extrapolation