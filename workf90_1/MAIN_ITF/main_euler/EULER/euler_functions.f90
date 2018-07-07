module euler_functions

! This module provides all the functions used in the Euler system.

implicit none

real(8), parameter :: error_gamma=-1.0d4

public  flux_in_x, flux_in_y, alarm_in_density, triem, uf, vf, &
        pf, cf, hf, inteng, tran3, tran5, alarm_in_pressure, roe, &
		state_error_message


contains


! The flux function in x-direction. The array w is the four 
! conserved variables, i.e. the density, the momentum 
! along the x-direction, the momentum along the y-direction, and 
! the total energy. f(j,w) is the jth component of the flux.
function flux_in_x(j,w,g)	result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: g
integer, intent(in) :: j
real(8) :: c
if(dabs(w(1)).lt.1.0d-12) then
 print*, 'Something must be wrong in density!'; pause
end if
select case(j)
 case(1); c=w(2)
 case(2); c=w(2)**2.0d0/w(1)+(g-1.0d0)* &
            (w(4)-0.5d0*(w(2)**2.0d0+w(3)**2.0d0)/w(1))
 case(3); c=w(2)*w(3)/w(1)
 case(4); c=w(2)/w(1)*(w(4)+(g-1.0d0)* &
            (w(4)-0.5d0*(w(2)**2.0d0+w(3)**2.0d0)/w(1)))
 case default; print*, 'Something must be wrong!!!'; pause
end	select
end function flux_in_x


! The flux function in y-direction. The array w is the same as in 
! the previous function f.
function flux_in_y(j,w,g)	result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: g
integer, intent(in) :: j
real(8) :: c
if(dabs(w(1)).lt.1.0d-12) then
 print*, 'Something must be wrong in density!'; pause
end if
select case(j)
 case(1); c=w(3)
 case(2); c=w(2)*w(3)/w(1)
 case(3); c=w(3)**2.0d0/w(1)+(g-1.0d0)* &
            (w(4)-0.5d0*(w(2)**2.0d0+w(3)**2.0d0)/w(1))
 case(4); c=w(3)/w(1)*(w(4)+(g-1.0d0)* &
            (w(4)-0.5d0*(w(2)**2.0d0+w(3)**2.0d0)/w(1)))
 case default; print*, 'Something must be wrong!!!'; pause
end	select
end function flux_in_y


! The function uf is the x component of the flow, which is 
! calculated by dividing the momentum along the X direction by 
! the density.
function uf(w) result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8) :: c
c=w(2)/w(1)
end	function uf


! The function vf is the y component of the flow, which is 
! calculated by dividing the momentum along the Y direction by 
! the density.
function vf(w) result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8) :: c
c=w(3)/w(1)
end function vf


! The function pf is the pressure in the flow, which is calculated 
! from the internal energy.
function pf(w,gamma) result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: gamma
real(8) :: c
if(gamma.lt.0.9d0*error_gamma) call state_error_message
c=(gamma-1.0d0)*(w(4)-0.5d0*(w(2)**2.0d0+w(3)**2.0d0)/w(1))
!if(c.lt.1.0d-12) then
! print*, 'Something must be wrong in PRESSURE!'; pause
!end if
end function pf


! The function cf is the sound speed of the flow, which is 
! calculated from the pressure and the density.
function cf(w,gamma) result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: gamma
real(8) :: c
c=dsqrt(gamma*pf(w,gamma)/w(1))
end function cf


! The function hf is the enthalpy of the flow, which is calculated 
! from the density, the pressure, and the total energy.
function hf(w,gamma) result(c)
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: gamma
real(8) :: c
c=(w(4)+pf(w,gamma))/w(1)
end function hf

! The function internal_energy is the internal energy of the flow, which is calculated 
! from the density, the pressure, and the total energy.
function inteng(w) result(c)
implicit none
real(8), dimension(4), intent(in) :: w
!real(8), intent(in) :: gamma
real(8) :: c
c=w(4)-0.5d0*(w(2)*w(2)+w(3)*w(3))/w(1)
end function inteng


! The function computes the pressure of a gamma-law fluid from its density and internal energy.
function pressure_fr_density_intengy(density,internal_energy,gamma) result(c)
implicit none
real(8), intent(in) :: density, internal_energy, gamma
real(8) :: c
if(gamma.lt.0.9d0*error_gamma) call state_error_message
c=(gamma-1.0d0)*internal_energy
end function pressure_fr_density_intengy


! Tran3 calculates a conserved array from a given row of density, 
! velocities and pressure.
subroutine tran3(rho,u,v,p,w,gamma)
implicit none
real(8), intent(in) :: rho, u, v, p
real(8), dimension(4), intent(out) :: w
real(8), intent(in) :: gamma
if(gamma.lt.0.9d0*error_gamma) call state_error_message
w(1)=rho
w(2)=u*rho
w(3)=v*rho
w(4)=p/(gamma-1.d0)+0.5d0*rho*(u**2.d0+v**2.d0)
end subroutine tran3


! tran5 calculates density, velocities, and pressure from a given 
! of conserved array.
subroutine tran5(rho,u,v,p,w,gamma)
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: gamma
real(8), intent(out) :: rho, u, v, p
rho=w(1)
u=uf(w)
v=vf(w)
p=pf(w,gamma)
end subroutine tran5


subroutine alarm_in_density(w)
! The subroutine 'alarm_in_density' checks the density and gives
! warning if the density is close to zero or negative.
implicit none
real(8), dimension(4), intent(in) :: w
if(w(1).lt.1.0d-12) then
 print*, 'Something must be wrong in DENSITY!'; pause
end if
end subroutine alarm_in_density


subroutine alarm_in_pressure(w,gamma)
! The subroutine 'alarm_in_pressure' checks the pressure and gives
! warning if the pressure is close to zero or negative.
implicit none
real(8), dimension(4), intent(in) :: w
real(8), intent(in) :: gamma
real(8) :: psu
psu=pf(w,gamma)
if(psu.lt.1.0d-12) then
 print*, 'Something must be wrong in PRESSURE!'; pause
end if
end subroutine alarm_in_pressure


! The subroutine triem solves the two dimensional Riemann problem.
subroutine triem(wl,wr,wml,wmr,gamma,al,be)
implicit none
real(8), dimension(4), intent(in) :: wl,wr
real(8), dimension(4), intent(out) :: wml,wmr
real(8), intent(in) :: gamma
real(8), intent(in) :: al, be

real(8), dimension(3) :: tl,tr
real(8)	:: rhol,ul, vl, pl, rhor, ur, vr, pr
real(8) :: uul, vvl, uur, vvr, qql, qqr, chl, chr
real(8) :: uml, vml, umr, vmr

real(8) :: gammus, beta, twa, twas

gammus=gamma-1.0d0
beta=(gamma+1.0d0)/gammus
twa=0.5d0*gammus/gamma
twas=dsqrt(twa)

call tran5(rhol,ul,vl,pl,wl,gamma)
call tran5(rhor,ur,vr,pr,wr,gamma)
uul=ul*al+vl*be
vvl=vl*al-ul*be
uur=ur*al+vr*be
vvr=vr*al-ur*be

qql=uul*uul+vvl*vvl
qqr=uur*uur+vvr*vvr
chl=0.5d0*rhol*uul*qql+rhol*uul*wl(4)+uul*pl
chr=0.5d0*rhor*uur*qqr+rhor*uur*wr(4)+uur*pr

call rsolv(rhol,uul,pl,rhor,uur,pr,tl,tr,gamma)
uml=tl(2)*al-vvl*be
vml=vvl*al+tl(2)*be
umr=tr(2)*al-vvr*be
vmr=vvr*al+tr(2)*be
call tran3(tl(1),uml,vml,tl(3),wml,gamma)
call tran3(tr(1),umr,vmr,tr(3),wmr,gamma)


contains


subroutine rsolv(dl,ul,pl,dr,ur,pr,tl,tr,gamma)
implicit none
real(8), intent(in) :: dl, ul, pl, dr, ur, pr, gamma
real(8), dimension(3), intent(out) :: tl, tr

real(8) :: cl, a, b, blog, bas, x1, xx1, xx2, xx3, yy1, dh1
real(8) :: dh2, c, zz1, zz2, x2, x3
integer :: i

cl=dsqrt(gamma*pl/dl)
a=dr/dl
b=pr/pl
blog=dlog(b)
c=(ur-ul)/cl
bas=dsqrt(b/a)
x1=0.
do i=1,100
 xx1=h1(x1)
 xx2=x1+blog
 xx3=h1(xx2)
 yy1=xx1+bas*xx3-c
 if (x1.ge.0.0d0) then
  dh1=2.0d0*twa*dexp(-twa*x1)/(gamma-1.0d0)
 else
  zz1=dexp(-x1)
  zz2=dsqrt(1.0d0+beta*zz1)
  dh1=2.0d0*twas/(gamma-1.0d0)*(zz1*zz2+(1.0d0-zz1)*	 &
  0.5d0*beta/zz2 )/zz2/zz2
 endif
 if (xx2.ge.0.0d0) then
  dh2=2.0d0*twa*dexp(-twa*xx2)/(gamma-1.0d0)
 else
  zz1=dexp(-xx2)
  zz2=dsqrt(1+beta*zz1)
  dh2=2.0d0*twas/(gamma-1)*(zz1*zz2+(1.0d0-zz1)*	 &
  0.5d0*beta/zz2 )/zz2/zz2
 endif
 x1=x1-yy1/(dh1+bas*dh2)
end do
x3=x1+blog
x2=dlog( a/f1(x1)*f1(x3) )
tl(1)=f1(x1)*dl
tl(2)=ul+cl*h1(x1)
tl(3)=dexp(-x1)*pl
tr(1)=dexp(x2)*tl(1)
tr(2)=tl(2)
tr(3)=tl(3)
!	t3d=tr(1)/f1(x3)
!	t3p=dexp(x3)*tr(2)
!	t3u=tr(3)+dsqrt(gamma*tr(2)/tr(1))*h3(x3)
end subroutine rsolv


function f1(x)
implicit none 
real(8), intent(in) :: x
real(8) :: f1
real(8) :: xx
if (x.ge.0.0d0) then
 f1=dexp(-x/gamma)
else
 xx=dexp(x)
 f1=(beta+xx)/(1.0d0+beta*xx)
endif
return
end	function f1


function h1(x)
implicit none
real(8), intent(in) :: x
real(8) :: h1
real(8) :: xx
if (x.ge.0.0d0) then
 h1=2.0d0*(1.0d0-dexp(-twa*x))/(gamma-1.0d0)
else
 xx=dexp(-x)
 h1=2.0d0*dsqrt(twa)*(1.0d0-xx)/(gamma-1.0d0)/dsqrt(1.0d0+beta*xx)
endif
return
end	function h1


function h3(x)
implicit none
real(8), intent(in) :: x
real(8) :: h3
real(8) :: xx
if (x.ge.0.0d0) then
 h3=2.0d0*(dexp(twa*x)-1.0d0)/(gamma-1.0d0)
else
 xx=dexp(x)
 h3=2.0d0*dsqrt(twa)*(xx-1.0d0)/(gamma-1.0d0)/dsqrt(1.0d0+beta*xx)
endif
return
end function h3


end subroutine triem


subroutine roe(wr,wl,el,er,ch,coef,g)

! Given two states w1, and w2 in density, momentums, and total 
! energy, the subroutine roe constructs an linearly approximate 
! matrix in the way suggested by Roe (cf. Roe's (1981)). The output 
! of this subroutine is the approximate eigenvalues char, the 
! approximate eigenvectors e, and the coefficients in the expression 
! of w1-w2 in e coef.

implicit none

real(8), dimension(4), intent(in) :: wl, wr
real(8), intent(in) :: g
real(8), dimension(4,4), intent(out) :: el, er
real(8), dimension(4), intent(out) :: ch, coef

real(8) :: rhol, ul, vl, pl, hl, rhor, ur, vr, pr, hr
real(8) :: rrho, uu, vv, hh, qq, aa, check

real(8), dimension(4) :: vec1, vec2
integer :: i

call tran5(rhol,ul,vl,pl,wl,g)
call tran5(rhor,ur,vr,pr,wr,g)
hl=hf(wl,g)
hr=hf(wr,g)
rrho=dsqrt(rhol)+dsqrt(rhor)
uu=(dsqrt(rhol)*ul+dsqrt(rhor)*ur)/rrho
vv=(dsqrt(rhol)*vl+dsqrt(rhor)*vr)/rrho
hh=(dsqrt(rhol)*hl+dsqrt(rhor)*hr)/rrho
qq=dsqrt(uu**2.0d0+vv**2.0d0)
aa=dsqrt((g-1.0d0)*(hh-0.5d0*qq**2.0d0))

! The left eigenvectors.
el(1,1)=uu*aa/(g-1.0d0)+0.5d0*(uu*uu+vv*vv)
el(1,2)=-uu-aa/(g-1.0d0)
el(1,3)=-vv
el(1,4)=1.0d0

el(2,1)=-vv
el(2,2)=0.0d0
el(2,3)=1.0d0
el(2,4)=0.0d0

el(3,1)=uu*uu-hh
el(3,2)=-uu
el(3,3)=0.0d0
el(3,4)=1.0d0

el(4,1)=-uu*aa/(g-1.0d0)+0.5d0*(uu*uu+vv*vv)
el(4,2)=-uu+aa/(g-1.0d0)
el(4,3)=-vv
el(4,4)=1.0d0

! The right eigenvectors.
er(1,1)=1.0d0
er(1,2)=uu-aa
er(1,3)=vv
er(1,4)=hh-uu*aa

er(2,1)=0.0d0
er(2,2)=0.0d0
er(2,3)=1.0d0
er(2,4)=vv

er(3,1)=1.0d0
er(3,2)=uu
er(3,3)=vv
er(3,4)=0.5d0*qq**2.0d0

er(4,1)=1.0d0
er(4,2)=uu+aa
er(4,3)=vv
er(4,4)=hh+uu*aa

! The eigenvalues.
ch(1)=uu-aa
ch(2)=uu
ch(3)=uu
ch(4)=uu+aa

coef(3)=(hh-qq**2.0d0)*(wr(1)-wl(1))+uu*(wr(2)-wl(2))+vv*(wr(3)-wl(3))
coef(3)=(coef(3)-(wr(4)-wl(4)))*(g-1.0d0)/aa**2.0d0
coef(2)=(wr(3)-wl(3))-vv*(wr(1)-wl(1))
coef(1)=0.5d0*((1.0d0+uu/aa)*(wr(1)-wl(1))-coef(3)-(wr(2)-wl(2))/aa)
coef(4)=0.5d0*((1.0d0-uu/aa)*(wr(1)-wl(1))-coef(3)+(wr(2)-wl(2))/aa)

vec1=el(1,:)
do i=2,4
 vec2=er(i,:)
 check=inner_product(vec1,vec2)
 if(dabs(check).gt.1.0d-8) then
  print*, 'There must be something wrong in Riemann approximate solver!!!'
  pause
 end if
end do

vec1=el(4,:)
do i=1,3
 vec2=er(i,:)
 check=inner_product(vec1,vec2)
 if(dabs(check).gt.1.0d-8) then
  print*, 'There must be something wrong in Riemann approximate solver!!!'
  pause
 end if
end do

do i=2,3
 vec1=el(i,:)
 vec2=er(1,:)
 check=inner_product(vec1,vec2)
 if(dabs(check).gt.1.0d-8) then
  print*, 'There must be something wrong in Riemann approximate solver!!!'
  pause
 end if
 vec2=er(4,:)
 check=inner_product(vec1,vec2)
 if(dabs(check).gt.1.0d-8) then
  print*, 'There must be something wrong in Riemann approximate solver!!!'
  pause
 end if
end do


contains


function inner_product(v1,v2) result(c)

implicit none
real(8), dimension(4), intent(in) :: v1, v2
real(8) :: c

integer :: i

c=0.0d0
do i=1,4
 c=c+v1(i)*v2(i)
end do

end function inner_product


end	subroutine roe


subroutine state_error_message

print*, 'Something is wrong here in physics!!!'
pause

end subroutine state_error_message



end module euler_functions