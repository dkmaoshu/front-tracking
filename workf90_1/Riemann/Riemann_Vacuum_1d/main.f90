Program Main

use Riemann_solver

implicit none

real*8::pl, pr, rhol, rhor, vl, vr
real*8::pl_star, rhol_star, vl_star
real*8::pr_star, rhor_star, vr_star
character*11,dimension(3) :: wv_pt

type(phys)::ur, ul, ur_star, ul_star, ull, urr

 type(phys) :: fr_star, fr
 real(8) :: sksp, ent, rrr, sss

ull%type_of_fluid='Polytropic          '
nullify(ull%vdw); allocate(ull%pol)
urr%type_of_fluid='Polytropic          '
nullify(urr%vdw); allocate(urr%pol)

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.125d0; vr=0.0d0; pr=0.1d0; urr%pol%thermo%gamma=1.4d0

rhor=1.0d0; vr=5.0d0; pr=1.0d0; urr%pol%thermo%gamma=1.4d0
rhol=1.0d0; vl=-5.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=1.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=0.125d0; vr=0.0d0; pr=1.0d0; br=0.0d0; gammr=1.4d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=1.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; br=0.0d0; gammr=1.4d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=0.125d0; vr=0.0d0; pr=0.1d0; br=0.0d0; gammr=1.2d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; gamml=7.0d0; bl=3000.0d0; rho0l=1.0d0
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

 ent=sf(ul)
 ent=sf(ul_star)
 rrr=vf(ul)+2.0d0*cf(ul)/0.4d0
 rrr=vf(ul_star)+2.0d0*cf(ul_star)/0.4d0

 ent=sf(ur)
 ent=sf(ur_star)
 sss=vf(ur)-2.0d0*cf(ur)/0.4d0
 sss=vf(ur_star)-2.0d0*cf(ur_star)/0.4d0

write(*,*)

write(*,*) 'The left meiddle state is'
write(*,'(3f25.16)') rhol_star,vl_star,pl_star
write(*,*) 'The right meiddle state is'
write(*,'(3f25.16)') rhor_star,vr_star,pr_star

End program Main
 
