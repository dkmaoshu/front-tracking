module produce_poly_nomial_4_xx_yy
! This module is for producing polynomial $p(x,y)$ on the left and right sides of a critical
! cell.

use interpolation
! 'interpl.f90'

use solu_comput
! 'solu_com.f90'

implicit none


public  p_side


contains


function p_side(temp,x,y,side) result(c)
! The polynomial function associated with the critical cell pointed by 'temp'. Symbolics in the
! remarks here can be refered to 'reset.tex'

implicit none

type(auxi_crit_cell), pointer :: temp
real(8), intent(in) :: x, y
character*6, intent(in) :: side
type(state) :: c

type(auxi_crit_cell), pointer :: temp_nb
integer, dimension(-1:1) :: i0, j0
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(state), dimension(-3:3) :: uxy
! The data sector.
character*3 :: xy_dir
! The direction of the preparing.
integer :: i
type(state), dimension(-1:1) :: uuxy
real(8), dimension(-1:1) :: xxyy
real(8) :: xy

g_cell=temp%g_cell; p_cell=temp%p_cell
if(g_cell%g_type.ne.'xx'.or.g_cell%g_type.ne.'yy') call error_message

! Produce reconstruction polynomial $p^{n,\pm}_j(x)$.
xxyy(-1)=-1.0d0; xxyy(0)=0.0d0; xxyy(1)=1.0d0
xy_dir=g_cell%g_type
i0(0)=g_cell%x_idx; j0(0)=g_cell%y_idx
call smth_dpr(i0(0),j0(0),temp,uxy,xy_dir,side,'critical',3)
select case(xy_dir)
 case('xx '); xy=x
 case('yy '); xy=y
end select
uuxy(0)=state_inter(uxy,xxyy,y)

! Produce reconstruction polynomials $p^{n,\pm}_{j+1}(x)$ and $p^{n,\pm}_{j-1}(x)$.
temp_nb=>temp%p_nxt_dr
g_cell=temp_nb%g_cell; p_cell=temp_nb%p_cell
i0(-1)=g_cell%x_idx; j0(-1)=g_cell%y_idx
call smth_dpr(i0(-1),j0(-1),temp,uxy,xy_dir,side,'critical',3)
select case(xy_dir)
 case('xx ') 
  if(i0(-1).gt.i0(0)) then
   xy=x-0.5d0
  else
   if(i0(-1).lt.i0(0)) then
    xy=x+0.5d0
   else
    xy=x
   end if
  end if
 case('yy ')
  if(j0(-1).gt.j0(0)) then
   xy=y-0.5d0
  else
   if(j0(-1).lt.j0(0)) then
    xy=y+0.5d0
   else
    xy=y
   end if
  end if
end select
uuxy(-1)=state_inter(uxy,xxyy,xy)

temp_nb=>temp%n_nxt_dr
g_cell=temp_nb%g_cell; p_cell=temp_nb%p_cell
i0(1)=g_cell%x_idx; j0(1)=g_cell%y_idx
call smth_dpr(i0(1),j0(1),temp,uxy,xy_dir,side,'critical',3)
select case(xy_dir)
 case('xx ') 
  if(i0(1).gt.i0(0)) then
   xy=x-0.5d0
  else
   if(i0(1).lt.i0(0)) then
    xy=x+0.5d0
   else
    xy=x
   end if
  end if
 case('yy ')
  if(j0(1).gt.j0(0)) then
   xy=y-0.5d0
  else
   if(j0(1).lt.j0(0)) then
    xy=y+0.5d0
   else
    xy=y
   end if
  end if
end select
uuxy(1)=state_inter(uxy,xxyy,xy)

! Produce $p^{n,-}(x,y)
select case(xy_dir)
 case('xx ')
  do i=-1,1
   xxyy(i)=dfloat(i0(i)-i0(0))
  end do
 case('yy ')
  do i=-1,1
   xxyy(i)=dfloat(j0(i)-j0(0))
  end do
end select
select case(xy_dir)
 case('xx '); xy=y
 case('yy '); xy=x
end select
c=state_inter(uuxy,xxyy,xy)


contains


function state_inter(uxy,xxyy,xy) result(c)
! The interpolation function via prime function on physical states.

implicit none

type(state), dimension(-1:1), intent(in) :: uxy
real(8), dimension(-1:1), intent(in) :: xxyy
real(8), intent(in) :: xy
type(state) :: c

integer :: i, j
real(8), dimension(3) :: vv, xyp
real(8) :: cc

do j=1,4
 do i=-1,1
  vv(i+1)=uxy(i)%value(j); xyp(i+1)=xxyy(i)
 end do
 cc=interpol_via_prime_fun(vv,xxyy,xy,3,0.5d0)
 c%value(j)=cc
end do

end function state_inter


end function p_side


end module produce_poly_nomial_4_xx_yy



