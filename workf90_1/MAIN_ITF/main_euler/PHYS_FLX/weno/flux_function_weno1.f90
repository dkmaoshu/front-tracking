module flux_functions
! This module provides all the functions used in the construction of
! the second-odrer TVD numerical flux function.

use grid
use euler_functions

implicit none

real(8), dimension(-2:3,4) :: f, wuu, fff
real(8), dimension(4,4) :: evl, evr, evll, evrr
real(8), dimension(4) :: emxy
!real(8), parameter :: gm1=gamma-1.0d0
real(8), parameter :: epweno=1.0d-6

public  wflux_x, wflux_y
private evl, evr, weno, evll, evrr, f,wuu, fff, epweno

contains


subroutine wflux_x(wu,flux,gamma)

! The subroutine wflux_x constructs a WENO flux in the x direction in the 
! way described in Guang-Shan Jiang and Chi-Wang Shu in Journal of 
! Computational Physics, v126 (1996), pp.202-228.

implicit none

real(8), dimension(-2:3,4), intent(in) :: wu
real(8), dimension(4), intent(out) :: flux 
real(8), intent(in) :: gamma

real(8), dimension(-2:3) :: w, vx, vy, h
integer :: i, j, kx
real(8) :: den, xmt, ymt, t0, eng, vex, vey, pre, cvel, t1, vxm, vym, hm, &
           qm, cm
real(8) :: rcm, b1, b2, t2, t3
real(8) :: gm1

 real(8) :: hhhh

gm1=gamma-1.0d0

do i=-2, 3
 do j=1, 4
  f(i,j)=flux_in_x(j,wu(i,:),gamma)
 end do
 den=wu(i,1)
 t0=1.0d0/den
 eng=wu(i,4)
 pre=pf(wu(i,:),gamma)
 w(i)=dsqrt(den)
 h(i)=( pre + eng ) * t0
 vx(i)=uf(wu(i,:))
 vy(i)=vf(wu(i,:))
end do

! Compute left and right eigenvectors of Roe's mean matrix.

t0 = w(0) / ( w(0) + w(1) )
t1 = 1.0d0 - t0
vxm = t0 * vx(0) + t1 * vx(1)
vym = t0 * vy(0) + t1 * vy(1)
hm = t0 *  h(0) + t1 *  h(1)
qm = 0.5d0 * ( vxm * vxm + vym * vym )
cm = dsqrt( gm1 * ( hm - qm ) )
t0 = vxm * cm
evr(1,1) = 1.0d0
evr(1,2) = 0.0d0
evr(1,3) = 1.0d0
evr(1,4) = 1.0d0
evr(2,1) = vxm - cm
evr(2,2) = 0.0d0
evr(2,3) = vxm
evr(2,4) = vxm + cm
evr(3,1) = vym
evr(3,2) = 1.0d0
evr(3,3) = vym
evr(3,4) = vym
evr(4,1) = hm - t0
evr(4,2) = vym
evr(4,3) = qm
evr(4,4) = hm + t0
rcm = 1.d0 / cm
b1 = gm1 * rcm * rcm
b2 = qm * b1
t0 = vxm * rcm
t1 = b1 * vxm
t2 = 0.5d0 * b1
t3 = b1 * vym
evl(1,1) = 0.5d0 * ( b2 + t0 )
evl(1,2) = -0.5d0 * ( t1 + rcm )
evl(1,3) = -0.5d0 * t3
evl(1,4) = t2
evl(2,1) = - vym
evl(2,2) = 0.0d0
evl(2,3) = 1.0d0
evl(2,4) = 0.0d0
evl(3,1) = 1.d0 - b2
evl(3,2) = t1
evl(3,3) = t3
evl(3,4) = -b1
evl(4,1) =  0.5d0 * ( b2 - t0 )
evl(4,2) = -0.5d0 * ( t1 - rcm )
evl(4,3) = -0.5d0 *   t3
evl(4,4) = t2

do i=-2, 3
 do j=1,4
  if(j.eq.2) then
   wuu(i,j)=-wu(-i+1,j)
  else
   wuu(i,j)=wu(-i+1,j)
  end if
 end do
end do
do i=-2, 3
 do j=1, 4
  fff(i,j)=flux_in_x(j,wuu(i,:),gamma)
 end do
end do

evll=evl; evll(2,:)=-evl(2,:); evll(1,:)=evl(4,:); evll(4,:)=evl(1,:)

evrr=evr; evrr(:,2)=-evr(:,2); evrr(:,1)=evr(:,4); evrr(:,4)=evr(:,1)

do i=1, 4
 do j=1, 4
  hhhh=0.0d0
  do kx=1, 4 
   hhhh=hhhh+evll(i,kx)*evrr(kx,j)
  end do
  continue
 end do
end do

! Call the WENO slover.

call weno(wu,2,flux)

return

end subroutine wflux_x


subroutine wflux_y(wu,flux,gamma)

! The subroutine wflux_x constructs a WENO flux in the y direction in the 
! way described in Guang-Shan Jiang and Chi-Wang Shu in Journal of 
! Computational Physics, v126 (1996), pp.202-228.

implicit none

real(8), dimension(-2:3,4), intent(in) :: wu
real(8), dimension(4), intent(out) :: flux 
real(8), intent(in) :: gamma

real(8), dimension(-2:3) :: w, vx, vy, h
integer :: i, j
real(8) :: den, xmt, ymt, t0, eng, vex, vey, pre, cvel, t1, vxm, vym, hm, &
           qm, cm
real(8) :: rcm, b1, b2, t2, t3
real(8) :: gm1

gm1=gamma-1.0d0

do i=-2, 3
 do j=1, 4
  f(i,j)=flux_in_y(j,wu(i,:),gamma)
 end do 
 den=wu(i,1)
 t0=1.0d0/den
 eng=wu(i,4)
 pre=pf(wu(i,:),gamma)
 w(i)=dsqrt(den)
 h(i)=( pre + eng ) * t0
 vx(i)=uf(wu(i,:))
 vy(i)=vf(wu(i,:)) 
end do

! Compute left and right eigenvectors of Roe's mean matrix.

t0 = w(0) / ( w(0) + w(1) )
t1 = 1.0d0 - t0
vxm = t0 * vx(0) + t1 * vx(1)
vym = t0 * vy(0) + t1 * vy(1)
hm = t0 *  h(0) + t1 *  h(1)
qm = 0.5d0 * ( vxm * vxm + vym * vym )
cm = dsqrt( gm1 * ( hm - qm ) )
t0 = vym * cm
evr(1,1) = 1.0d0
evr(1,2) = 0.0d0
evr(1,3) = 1.0d0
evr(1,4) = 1.0d0
evr(2,1) = vxm
evr(2,2) = 1.0d0
evr(2,3) = vxm
evr(2,4) = vxm
evr(3,1) = vym - cm
evr(3,2) = 0.0d0
evr(3,3) = vym
evr(3,4) = vym + cm
evr(4,1) = hm - t0
evr(4,2) = vxm
evr(4,3) = qm
evr(4,4) = hm + t0
rcm = 1.0d0 / cm
b1 = gm1 * rcm * rcm
b2 = qm * b1
t0 = vym * rcm
t1 = b1 * vym
t2 = 0.5d0 * b1
t3 = b1 * vxm
evl(1,1) = 0.5d0 * ( b2 + t0 )
evl(1,2) = -0.5d0 * t3
evl(1,3) = -0.5d0 * ( t1 + rcm )
evl(1,4) = t2
evl(2,1) = - vxm
evl(2,2) = 1.0d0
evl(2,3) = 0.0d0
evl(2,4) = 0.0d0
evl(3,1) = 1.0d0 - b2
evl(3,2) = t3
evl(3,3) = t1
evl(3,4) = -b1
evl(4,1) =  0.5d0 * ( b2 - t0 )
evl(4,2) = -0.5d0 *   t3
evl(4,3) = -0.5d0 * ( t1 - rcm )
evl(4,4) = t2

do i=-2, 3
 do j=1,4
  if(j.eq.3) then
   wuu(i,j)=-wu(-i+1,j)
  else
   wuu(i,j)=wu(-i+1,j)
  end if
 end do
end do
do i=-2, 3
 do j=1, 4
  fff(i,j)=flux_in_y(j,wuu(i,:),gamma)
 end do
end do

evll=evl; evll(3,:)=-evl(3,:)
evrr=evr; evrr(:,3)=-evr(:,3)

! Call the WENO solver.

call weno(wu,3,flux)

return

end subroutine wflux_y


subroutine weno(wu,xy_dir,flux)

implicit none

!****678****************************************************************       
! Name:	     wenolf.f
! Function:  Use WENO-LF-4 or WENO-LF-5 to approximate fluxes 
!****678****************************************************************       

real(8), dimension(-2:3,4), intent(in) :: wu
integer, intent(in) :: xy_dir
real(8), dimension(4), intent(out) :: flux

real(8), dimension(-2:2,4) :: df, du, dff, duu
real(8), dimension(4) :: ff
real(8), dimension(-2:2,4,2) :: gg
real(8), dimension(4,2) :: hh
integer :: m, i, nsm, m1, k0, k1
real(8) :: t1, t2, t3, tt1, tt2, tt3, s1, s2, s3, t0, llll, rrrr

!real(8), dimension(4,4) :: evll
!real(8), dimension(4,2) :: hhh
real(8), dimension(4,2) :: zzzz

! real(8), dimension(-2:2,4) :: ggxx
! real(8) :: xxxx, yyyy

!nsm = ns + md

!evll=evl; evll(:,2)=-evll(:,2)

do m=1,4
 do i=-2, 2
  df(i,m)=f(i+1,m)-f(i,m)
  du(i,m)=wu(i+1,m)-wu(i,m)
  dff(i,m)=fff(i+1,m)-fff(i,m)
  duu(i,m)=wuu(i+1,m)-wuu(i,m)
 enddo
enddo

!----------------- Loop in "m" starts here.  -------------------

do m=1, 4

! Use Lax-Friedrichs building block to split the fluxes

 do m1=1, 4
  do  i = -2, 2
   gg(i,m1,1) = 0.5d0 * ( df(i,m1) + emxy(m) * du(i,m1) )
!   gg(i,m1,2) = gg(i,m1,1) - df(i,m1)
   gg(i,m1,2) = 0.5d0 * ( dff(i,m1) + emxy(m) * duu(i,m1) )
  enddo
 enddo

!  ggxx=gg(:,:,2)

! Project the positive and negative part of the fluxes to the 'm'th 
! characteristic field

! do m1 = 1, 4
!  k0 = m1 - 3
!  k1 =  3 - m1
!  hh(m1,1)= evl(m,1)*gg(k0,1,1) + evl(m,2)*gg(k0,2,1) &
!            + evl(m,3)*gg(k0,3,1) + evl(m,4)*gg(k0,4,1)
!  hh(m1,2)= evl(m,1)*gg(k1,1,2) + evl(m,2)*gg(k1,2,2) &
!            + evl(m,3)*gg(k1,3,2) + evl(m,4)*gg(k1,4,2)
!           
          
! enddo

 do k0 = -2, 1
  hh(k0+3,1)= evl(m,1)*gg(k0,1,1) + evl(m,2)*gg(k0,2,1) + evl(m,3)*gg(k0,3,1) + evl(m,4)*gg(k0,4,1)
  hh(k0+3,2)= evll(m,1)*gg(k0,1,2) + evll(m,2)*gg(k0,2,2) + evll(m,3)*gg(k0,3,2) + evll(m,4)*gg(k0,4,2)
 end do
     
! Compute the weights and approximate the fluxes.
   
 ff(m) = 0.0d0
  
 do m1 = 1, 2
 
  t1 = hh(1,m1) - hh(2,m1)
  t2 = hh(2,m1) - hh(3,m1)
  t3 = hh(3,m1) - hh(4,m1)
   
! un-comment(comment) the following 3 lines to use (not to use) the new 
! smoothness measurement:
   
  if(m1.eq.1) then
   tt1 = 13.0d0*t1**2.0d0+3.0d0*(hh(1,m1)-3.0d0*hh(2,m1))**2.0d0
   tt2 = 13.0d0*t2**2.0d0+3.0d0*(hh(2,m1)+hh(3,m1))**2.0d0
   tt3 = 13.0d0*t3**2.0d0+3.0d0*(3.0d0*hh(3,m1)-hh(4,m1))**2.0d0
  else
   tt1 = 13.0d0*t1**2.0d0+3.0d0*(hh(1,m1)-3.0d0*hh(2,m1))**2.0d0
   tt2 = 13.0d0*t2**2.0d0+3.0d0*(hh(2,m1)+hh(3,m1))**2.0d0
   tt3 = 13.0d0*t3**2.0d0+3.0d0*(3.0d0*hh(3,m1)-hh(4,m1))**2.0d0
  end if
   
  tt1 =  ( epweno + tt1 )**2.0d0
  tt2 =  ( epweno + tt2 )**2.0d0
  tt3 =  ( epweno + tt3 )**2.0d0
  s1 =      tt2 * tt3
  s2 = 6.0d0 * tt1 * tt3
  s3 = 3.0d0 * tt1 * tt2
  t0 = 1.0d0 / ( s1 + s2 + s3 )
  s1 = s1 * t0
  s3 = s3 * t0
!  ff(m) = ff(m) + ( s1*(t2-t1) + (0.5d0*s3-0.25d0)*(t3-t2) ) /3.0d0
   
  zzzz(m,m1)=( s1*(t2-t1) + (0.5d0*s3-0.25d0)*(t3-t2) ) /3.0d0
    
 enddo
 
enddo
!----------------- Loop in "m"  ends  here.  -------------------

! Project the fluxes to the physical space:

do m =  1, 4
!  flux(m)= evr(m,1) * ff(1) + evr(m,2) * ff(2) &
!           + evr(m,3) * ff(3) + evr(m,4) * ff(4) &
!           + ( -f(-1,m) + 7.0d0*( f(0,m)+f(1,m) ) - f(2,m) )/12.0d0

 llll=evr(m,1)*zzzz(1,1)+evr(m,2)*zzzz(2,1)+evr(m,3)*zzzz(3,1)+evr(m,4)*zzzz(4,1)
 rrrr=evrr(m,1)*zzzz(1,2)+evrr(m,2)*zzzz(2,2)+evrr(m,3)*zzzz(3,2)+evrr(m,4)*zzzz(4,2)
  
 if(m.eq.xy_dir) then
  flux(m)= ( -f(-1,m) + 7.0d0*( f(0,m)+f(1,m) ) - f(2,m) )/12.0d0+llll+rrrr
 else
  flux(m)= ( -f(-1,m) + 7.0d0*( f(0,m)+f(1,m) ) - f(2,m) )/12.0d0+llll-rrrr
 end if

! xxxx=evr(m,1) * ff(1) + evr(m,2) * ff(2)+ evr(m,3) * ff(3) + evr(m,4) * ff(4)
! yyyy=-f(-1,m) + 7.0d0*( f(0,m)+f(1,m) ) - f(2,m) 
          
enddo

return

end subroutine weno


end module flux_functions