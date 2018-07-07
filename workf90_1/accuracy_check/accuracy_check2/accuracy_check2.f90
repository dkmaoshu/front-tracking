program accuracy_check2

implicit none

real(8) :: x_width, y_width, sum, sum1, xx, xi, yj, h
integer :: nxll, nyll, nx, ny
integer :: nxll_exact, nyll_exact, nx_exact, ny_exact
integer :: i, j, ii, jj, ic, jc, mod
!real(8) :: h, hc, errorf, errorc, order
real(8), dimension(:,:), allocatable :: u_show, u_show_1, u_exact

open(1,file='d:\workf90_1\accuracy_check\accuracy_check2\show\shown.dat')
read(1,'(e25.15)') x_width, y_width
read(1,'(i5)') nxll_exact, nyll_exact, nx_exact, ny_exact
close(1)

allocate(u_exact(nxll_exact:nxll_exact+nx_exact-1,nyll_exact:nyll_exact+ny_exact-1))

open(1, file='d:\workf90_1\accuracy_check\accuracy_check2\show\showrho.dat')
do i=nxll_exact,nxll_exact+nx_exact-1
 do j=nyll_exact,nyll_exact+ny_exact-1
  read(1,'(e25.15)') u_exact(i,j)
 end do
end do  
close(1)

open(1,file='d:\workf90_1\show\shown.dat')
read(1,'(e25.15)') x_width, y_width
read(1,'(i5)') nxll, nyll, nx, ny
close(1)

h=x_width/dfloat(nx-1)

allocate(u_show(nxll:nxll+nx-1,nyll:nyll+ny-1))

open(1, file='d:\workf90_1\show\showrho.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  read(1,'(e25.15)') u_show(i,j)
 end do
end do  
close(1)

mod=nxll_exact/nxll/2

sum=0.0d0
do i=nxll+1,nxll+nx-2
 do j=nyll+1,nyll+ny-2
  xi=i*h; yj=j*h
  sum1=0.0d0
  ic=2*mod*i; jc=2*mod*j
  do ii=-mod+1,mod-1
   do jj=-mod+1,mod-1
    sum1=sum1+u_exact(ic+ii,jc+jj)
   end do
  end do
  do ii=-mod+1,mod-1
   xx=0.5d0*(u_exact(ic+ii,jc-mod)+u_exact(ic+ii,jc+mod))
   sum1=sum1+xx
  end do
  do jj=-mod+1,mod-1
   xx=0.5d0*(u_exact(ic-mod,jc+jj)+u_exact(ic-mod,jc+jj))
   sum1=sum1+xx
  end do
  xx=0.25d0*(u_exact(ic-mod,jc-mod)+u_exact(ic+mod,jc-mod))
  xx=xx+0.25d0*(u_exact(ic-mod,jc+mod)+u_exact(ic+mod,jc+mod))
  sum1=sum1+xx
  sum1=sum1/dfloat(mod)/dfloat(mod)/4.0d0
  xx=dabs(u_show(i,j)-sum1)
  sum=sum+xx
 end do
end do

sum=sum/dfloat(nxll)/dfloat(nxll)/4.0d0

!xx=u_show(12,12)-u_exact(48,48)

!order=(dlog10(errorc)-dlog10(errorf))/dlog10(2.0d0)

continue 

end program accuracy_check2