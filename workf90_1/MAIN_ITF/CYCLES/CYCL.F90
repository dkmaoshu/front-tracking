program cycles

use grid
! 'grid.f90'

implicit none

real(8), dimension(0:100000) :: x, y
real(8) :: center_x, center_y, alpha, pi, xx, yy, position
integer :: i, steps, x_index, y_index
character*3 :: compute_position, direction
real(8), dimension(2) :: normal

print*, 'Input steps'
read(*,*) steps
print*, 'Input mesh-size'
read(*,*) h

pi=4.0d0*datan(1.0d0)
center_x=dfloat(steps)*r*h
center_y=dfloat(steps)*r*h
do while(center_x.gt.0.2d0) 
 center_x=center_x-0.2d0
 center_y=center_y-0.2d0
end do
do i=0,100000
 alpha=pi/50000.0d0*dfloat(i)
 x(i)=dcos(alpha)*0.6d0+center_x
 y(i)=dsin(alpha)*0.6d0+center_y
end do

open(1,file='d:\workf90_1\show\cycle.dat')
do i=0,100000
 write(1,'(2f12.6)') x(i), y(i)
end do
close(1)

print*, 'Do you want compute a discontinuity position?'
read(*,*) compute_position
do while(compute_position.eq.'yes') 
 print*, 'Please input the direction'
 read(*,*) direction
 select case(direction)
  case('xx ')
   print*, 'Please input the x-index'
   read(*,*) x_index
   print*, 'Please input the y-variable'
   read(*,*) yy
   position=center_x+dsqrt(0.36d0-(yy-center_y)*(yy-center_y))
   normal(1)=(position-center_x)/0.6d0
   normal(2)=(yy-center_y)/0.6d0
   position=position/h-dfloat(x_index)
  case('yy ')
   print*, 'Please input the y-index'
   read(*,*) y_index
   print*, 'Please input the x-variable'
   read(*,*) xx
   position=center_y+dsqrt(0.36d0-(xx-center_x)*(xx-center_x))
   normal(1)=(xx-center_x)/0.6d0
   normal(2)=(position-center_y)/0.6d0
   position=position/h-dfloat(y_index)
  case default; print*, 'The direction is wrong.'
 end select
 write(*,'(f16.6)') position
 write(*,'(2f16.6)') normal(1), normal(2)
 print*, 'Do you want compute one more discontinuity position?'
 read(*,*) compute_position
end do

end program cycles