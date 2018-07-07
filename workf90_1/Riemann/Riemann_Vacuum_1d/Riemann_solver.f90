module Riemann_solver

use fluids

implicit none


contains


subroutine riemann(ul,ur,ul_star,ur_star,wave_pattern)

implicit none

type(phys), intent(in) :: ul, ur
type(phys), intent(out) :: ul_star, ur_star
character*11, dimension(3) :: wave_pattern

real*8 :: ml, mr, ml_q, mr_q, p_q
real*8 :: p_star, v_star, rhol_star, rhor_star 
real*8 :: pl, pr, rhol, rhor, vl, vr
real*8 :: sl, sr
real*8 :: alf, epslon
real*8 :: al, ar
real*8 :: pp
real*8 :: difml, difmr
real*8 :: max
integer :: q, q_stop

 real(8) :: xxx

ml_q=error_data
mr_q=error_data
 
call tran5(rhol,vl,pl,ul) 
call tran5(rhor,vr,pr,ur)

epslon=1.0d-12
q=0
ml=100.0d0
mr=100.0d0

p_star=(pl+pr)*0.5d0

!coefl=-dsqrt((pl+bl)*rhol)
!coefr=dsqrt((pr+br)*rhor)
alf=1.0d0
p_q=0.d0
q_stop=100
   
ml_q=phi(rhol,pl,ul,p_star)
mr_q=phi(rhor,pr,ur,p_star)
difml=dabs(ml_q-ml)
difmr=dabs(mr_q-mr)

max=dmax1(difml,difmr)
pp=dabs(p_star-p_q)

do while(max.ge.epslon.and.pp.ge.epslon.and.q.lt.q_stop)
! if(p_star.lt.epslon)then
!  p_star=epslon
! end if
 ml=ml_q
 mr=mr_q
   
 p_q=p_star
 p_star=(vl-vr+pl/ml+pr/mr)/(1.0d0/ml+1.0d0/mr)
 p_star=alf*p_star+(1.0d0-alf)*p_q

 if(p_star.lt.epslon)then
  p_star=epslon
 end if
    
 ml_q=phi(rhol,pl,ul,p_star)
 mr_q=phi(rhor,pr,ur,p_star)
 difml=dabs(ml_q-ml)
 difmr=dabs(mr_q-mr)

  xxx=0.2/dsqrt(1.4d0)
   
 max=dmax1(difml,difmr)
 pp=dabs(p_star-p_q)
 q=q+1
 
end do

v_star=(pl-pr+ml*vl+mr*vr)/(ml+mr)
 
if(p_star.gt.pl)then
! The left center wave is a shock.
 sl=(rhol*vl-ml)/rhol
 rhol_star=ml/(v_star-sl)
 wave_pattern(1)='shock      '
 write(*,*)'The left center wave is a shock'
 else
!The left center wave is a rarefaction wave.
 rhol_star=density_rare(rhol,pl,ul,p_star)
 wave_pattern(1)='rarefaction'
 write(*,*)'The left center wave is a rarefaction wave'
end if

if(p_star.gt.pr)then
!The right center wave is a shock.
 sr=(rhor*vr+mr)/rhor
 rhor_star=-mr/(v_star-sr) 
 wave_pattern(3)='shock      '
 write(*,*)'The right center wave is a shock'
else
!The right center wave is a rarefaction wave.
 rhor_star=density_rare(rhor,pr,ur,p_star)
 write(*,*)'The right center wave is a rarefaction wave'
end if

wave_pattern(2)='contact   '

call phys_clon(ul_star,ul)
call phys_clon(ur_star,ur)

call tran3(rhol_star,v_star,p_star,ul_star)
call tran3(rhor_star,v_star,p_star,ur_star)
 
end subroutine riemann


end module Riemann_solver