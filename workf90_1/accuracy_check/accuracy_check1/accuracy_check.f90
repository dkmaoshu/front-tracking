program accuracy_check

implicit none

real(8) :: x_width, y_width
integer :: nxll, nyll, nx, ny
integer :: i, j
real(8) :: h, hc, errorf, errorc, order
real(8), dimension(:,:), allocatable :: u_show, u_show_1


open(1, file='d:\workf90_1\showf\shown.dat')
read(1,'(f10.4)') x_width, y_width
read(1,'(i5)') nxll, nyll, nx, ny
close(1)

allocate(u_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
open(2, file='d:\workf90_1\showf\showu.dat')
do i=nxll,nxll+nx-1; do j=nyll,nyll+ny-1
 read(2,'(f12.6)') u_show(i,j)
end do; end do  
close(2)

allocate(u_show_1(nxll:nxll+nx-1,nyll:nyll+ny-1))
open(2, file='d:\workf90_1\showf_1\showu.dat')
do i=nxll,nxll+nx-1; do j=nyll,nyll+ny-1
 read(2,'(f12.6)') u_show_1(i,j)
end do; end do  
close(2)

h=x_width/dfloat(nx-1)

errorf=0.0d0
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  errorf=errorf+dabs(u_show(i,j)-u_show_1(i,j))
 end do
end do
errorf=0.25d0*errorf*h*h

deallocate(u_show); deallocate(u_show_1)


open(1, file='d:\workf90_1\showc\shown.dat')
read(1,'(f10.4)') x_width, y_width
read(1,'(i5)') nxll, nyll, nx, ny
close(1)

allocate(u_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
open(2, file='d:\workf90_1\showc\showu.dat')
do i=nxll,nxll+nx-1; do j=nyll,nyll+ny-1
 read(2,'(f12.6)') u_show(i,j)
end do; end do  
close(2)

allocate(u_show_1(nxll:nxll+nx-1,nyll:nyll+ny-1))
open(2, file='d:\workf90_1\showc_1\showu.dat')
do i=nxll,nxll+nx-1; do j=nyll,nyll+ny-1
 read(2,'(f12.6)') u_show_1(i,j)
end do; end do  
close(2)

hc=x_width/dfloat(nx-1)

errorc=0.0d0
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  errorc=errorc+dabs(u_show(i,j)-u_show_1(i,j))
 end do
end do
errorc=0.25d0*errorc*hc*hc


deallocate(u_show); deallocate(u_show_1)

order=(dlog10(errorc)-dlog10(errorf))/dlog10(2.0d0)

continue 

end program accuracy_check