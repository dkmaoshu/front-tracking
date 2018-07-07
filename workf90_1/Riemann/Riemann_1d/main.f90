
Program Main

use Riemann_solver

implicit none

real*8::pl, pr, rhol, rhor, vl, vr
real*8::pl_star, rhol_star, vl_star
real*8::pr_star, rhor_star, vr_star
character*11,dimension(3) :: wv_pt

type(phys)::ur, ul, ur_star, ul_star, ull, urr, xxx

type(phys) :: fr_star, fr, fl, fl_star
real(8) :: sksp, ent

real(8) :: rhox, ux, px, sx, sxl, sxr, sr
real(8) :: aaa

!The kind of gas on the left.

ull%type_of_fluid='Polytropic          '
nullify(ull%vdw); nullify(ull%bar); allocate(ull%pol)

!ull%type_of_fluid='Barotropic          '
!nullify(ull%vdw); nullify(ull%pol); allocate(ull%bar)

!ull%type_of_fluid='Van Der Waals       '
!nullify(ull%pol); nullify(ull%bar); allocate(ull%vdw)

!The kind of gas on the right.

urr%type_of_fluid='Polytropic          '
nullify(urr%vdw); nullify(urr%bar); allocate(urr%pol)

!urr%type_of_fluid='Barotropic          '
!nullify(urr%vdw); nullify(urr%pol); allocate(urr%bar)

!urr%type_of_fluid='Van Der Waals       '
!nullify(urr%bar); nullify(urr%pol); allocate(urr%vdw)
  

!rhol=1.1d0; vl=0.0d0; pl=2000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0;
!ull%bar%thermo%p0=1.0d0
!rhor=1.1d0; vr=0.0d0; pr=1000.0d0; urr%bar%thermo%gamma=7.0d0; 
!urr%bar%thermo%beta=3000.0d0; urr%bar%thermo%rho0=1.0d0;
!urr%bar%thermo%p0=1.0d0

!rhol=5.0d1; vl=1.0d3; pl=1.0d+5; ull%vdw%thermo%gamma=1.4d0; 
!ull%vdw%thermo%a=5.0d0; ull%vdw%thermo%b=1.0d-3
!rhor=5.0d1; vr=0.0d0; pr=1.0d+5; urr%vdw%thermo%gamma=1.4d0; 
!urr%vdw%thermo%a=5.0d0; urr%vdw%thermo%b=1.0d-3

!rhol=2.27550882502480d0; vl=0.781011006283816d0; pl=2.56167433147329d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.167d0; vr=0.0d0; pr=1.013d0; urr%pol%thermo%gamma=1.63d0; 

!rhol=5.0d+1; vl=0.0d0; pl=1.0d+5; ull%vdw%thermo%gamma=1.4d0; 
!ull%vdw%thermo%a=5.0d0; ull%vdw%thermo%b=1.0d-3
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; urr%pol%thermo%gamma=1.4d0

!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0;
!ull%bar%thermo%p0=1.0d0
!rhor=5.0d+1; vr=0.0d0; pr=1.0d+5; urr%vdw%thermo%gamma=1.4d0; 
!urr%vdw%thermo%a=5.0d0; urr%vdw%thermo%b=1.0d-3

!rhol=5.0d+1; vl=0.0d0; pl=1.0d+5; ull%vdw%thermo%gamma=1.4d0; 
!ull%vdw%thermo%a=5.0d0; ull%vdw%thermo%b=1.0d-3
!rhor=1.0d0; vr=6.3386d0; pr=1000.0d0; urr%bar%thermo%gamma=7.0d0; 
!urr%bar%thermo%beta=3000.0d0; urr%bar%thermo%rho0=1.0d0;
!urr%bar%thermo%p0=1.0d0

!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0;
!ull%bar%thermo%p0=1.0d0
!rhor=0.125d0; vr=0.0d0; pr=0.1d0; urr%pol%thermo%gamma=1.2d0

!rhol=0.344568d0; vl=1.52872d0; pl=2.46610d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.5d0; vr=0.0d0; pr=0.571d0; urr%pol%thermo%gamma=1.4d0
!rr%bar%thermo%beta=3000.0d0; urr%bar%thermo%rho0=1.0d0;
!urr%bar%thermo%p0=1.0d0


! The Sod's shock tube problem for gamma-law gas.

!rhor=0.125d0; vr=0.0d0; pr=0.1d0; urr%pol%thermo%gamma=1.4d0
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0

!rhor=1.0d0; vr=0.0d0; pr=1.0d-2; urr%pol%thermo%gamma=1.4d0
!rhol=1.0d0; vl=0.0d0; pl=1.0d+3; ull%pol%thermo%gamma=1.4d0

!rhol=19.2374d0; vl=2000.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=1.29d-3; vr=0.0d0; pr=1.0d0; urr%pol%thermo%gamma=1.4d0

!rhol=115.0d0; vl=1000.0d0; pl=23084402d0; urr%pol%thermo%gamma=1.4d0
!rhor=1.290d-3; vr=0.0d0; pr=1.0d0; ull%pol%thermo%gamma=1.4d0

rhol=4.52d0; vl=0.0d0; pl=1.0d0; urr%pol%thermo%gamma=1.4d0
rhor=1.51569506726457d0; vr=-0.523345519274197d0; pr=1.805d0; ull%pol%thermo%gamma=1.4d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=1.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=0.125d0; vr=0.0d0; pr=1.0d0; br=0.0d0; gammr=1.4d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=1.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; br=0.0d0; gammr=1.4d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.125d0; vr=0.0d0; pr=0.1d0; urr%pol%thermo%gamma=1.4d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; gammr=7.0d0; br=3000.0d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=0.001d0; vl=0.0d0; pl=1000.0d0; gamml=1.4d0; bl=0.0d0; rho0l=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; gammr=7.0d0; br=3000.0d0; rho0r=1.0d0

! The Lax's shock tube problem for gamma-law gas.
!rhol=0.445d0; vl=0.698d0; pl=3.528d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.5d0; vr=0.0d0; pr=0.571d0; urr%pol%thermo%gamma=1.4d0

call phys_clon(ul,ull)
call tran3(rho=rhol,v=vl,p=pl,u=ul)
call phys_clon(ur,urr)
call tran3(rho=rhor,v=vr,p=pr,u=ur)

call riemann(ul,ur,ul_star,ur_star,wv_pt)

call tran5(rhol_star,vl_star,pl_star,ul_star)
call tran5(rhor_star,vr_star,pr_star,ur_star)

! aaa=-dlog(0.01d0)

 fl_star=f(ul_star) 
 fl=f(ul)

 sksp=(df(fl_star)-df(fl))/(df(ul_star)-df(ul))
 sksp=(mf(fl_star)-mf(fl))/(mf(ul_star)-mf(ul))
 sksp=(ef(fl_star)-ef(fl))/(ef(ul_star)-ef(ul))

 fl=f(ul)  
 fl_star=f(ul_star)

! xxx%type_of_fluid='Polytropic          '
! nullify(xxx%vdw); nullify(xxx%bar); allocate(xxx%pol)
! xxx%pol%thermo%gamma=1.4d0

! xxx%pol%csq%value(1)=ur%pol%csq%value(1)-0.2d0*(fr%pol%csq%value(1)-fl%pol%csq%value(1))
! xxx%pol%csq%value(2)=ur%pol%csq%value(1)-0.2d0*(fr%pol%csq%value(2)-fl%pol%csq%value(2))
! xxx%pol%csq%value(3)=ur%pol%csq%value(3)-0.2d0*(fr%pol%csq%value(3)-fl%pol%csq%value(3))

! rhox=df(xxx) 
! ux=vf(xxx) 
! px=pf(xxx)
! sx=sf(xxx)
 
! sxl=sf(ul_star); sxr=sf(ur_star); sr=sf(ur)

! sksp=(df(fl_star)-df(fr))/(df(ul)-df(ul_star))
! sksp=(mf(fl_star)-mf(fr))/(mf(ul)-mf(ul_star))
! sksp=(ef(fl_star)-ef(fr))/(ef(ul)-ef(ul_star))

! ent=sf(ul)
! ent=sf(ul_star)

! aaa=(px/dexp(2.272d0))**(1.0d0/1.4d0)
! aaa=(px/dexp(0.409d0))**(1.0d0/1.4d0)

! aaa=dlog(831.601d0)
! aaa=dlog(1.615d0)

write(*,*)

write(*,*) 'The left middle state is'
write(*,'(3f25.14)') rhol_star,vl_star,pl_star
write(*,*) 'The right middle state is'
write(*,'(3f25.14)') rhor_star,vr_star,pr_star

! ent=dlog(460.0d0)
! ent=(dlog(5.4d0)-dlog(4.4d0))*1.4d0
! ent=(dlog(4.4d0)-dlog(1.7d0))*1.4d0
 
End program Main
 
