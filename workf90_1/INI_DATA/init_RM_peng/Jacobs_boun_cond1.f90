module Jacobs_Peng_Zabusky_boundary

use solu_comput
! 'solu_comp.f90'

use numerical_fluxes
! 'flux.f90'

implicit none

public  Jacobs_PengZ_1, Jacobs_PengZ_2, Jacobs_PengZ_3, Jacobs_PZ_nonreflecting_1, Jacobs_PZ_nonreflecting_2

private smooth_region, smooth_rg_nonreflecting_1, smooth_rg_nonreflecting_2
 

contains


subroutine Jacobs_PengZ_1

implicit none

call smooth_region

call material_interface


contains


subroutine material_interface

implicit none

type(auxi_crit_cell), pointer :: temp_l, temp_r, current, temp_b
integer :: i, j, x_idx, y_idx, x_ii, y_jj
logical :: done, head_left

if(acvv(1)%status.eq.'asleep') return

! Determine whether the head of the interface is on the left or right.
i=acvv(1)%begin%next%g_cell%x_idx
if(i.eq.nxll) then
 head_left=.true.
else
 if(i.eq.nxll+nx-1) then
  head_left=.false.
 else
  call error_message
 end if
end if

! Periodic extension at the left end.
if(head_left) then
 temp_l=>acvv(1)%begin%next; temp_r=>acvv(1)%eend%previous%previous
else
 temp_r=>acvv(1)%begin%next%next; temp_l=>acvv(1)%eend%previous
end if

done=.false.
do while(.not.done)
 allocate(current)
 current=temp_r
 current%g_cell%x_idx=current%g_cell%x_idx-nx+1
 if(head_left) then
  call insert_boundary_cell(temp_l,current,'before')
 else
  call insert_boundary_cell(temp_l,current,'after ')
 end if
 x_idx=current%g_cell%x_idx; y_idx=current%g_cell%y_idx
 ggd_cell(x_idx,y_idx)%region='crit  '
 x_ii=temp_r%g_cell%x_idx; y_jj=temp_r%g_cell%y_idx
 do j=1, 4
  ggd_cell(x_idx,y_idx)%ccpt(j)%address=current%address
  ggd_cell(x_idx,y_idx)%ccpt(j)%side=ggd_cell(x_ii,y_jj)%ccpt(j)%side
 end do
 if(head_left) then
  temp_r=>temp_r%previous; temp_l=>temp_l%previous
 else
  temp_r=>temp_r%next; temp_l=>temp_l%next
 end if
 x_idx=temp_r%g_cell%x_idx
 if(x_idx.lt.nxll+nx-4) done=.true.
end do

! Periodic extension at the right end.
if(head_left) then
 temp_l=>acvv(1)%begin%next%next; temp_r=>acvv(1)%eend%previous
else
 temp_r=>acvv(1)%begin%next; temp_l=>acvv(1)%eend%previous%previous
end if

done=.false.
do while(.not.done)
 allocate(current)
 current=temp_l
 current%g_cell%x_idx=current%g_cell%x_idx+nx-1
 if(head_left) then
  call insert_boundary_cell(temp_r,current,'after ')
 else
  call insert_boundary_cell(temp_r,current,'before')
 end if
 x_idx=current%g_cell%x_idx; y_idx=current%g_cell%y_idx
 ggd_cell(x_idx,y_idx)%region='crit  '
 x_ii=temp_l%g_cell%x_idx; y_jj=temp_l%g_cell%y_idx
 do j=1, 4
  ggd_cell(x_idx,y_idx)%ccpt(j)%address=current%address
  ggd_cell(x_idx,y_idx)%ccpt(j)%side=ggd_cell(x_ii,y_jj)%ccpt(j)%side
 end do
 if(head_left) then
  temp_r=>temp_r%next; temp_l=>temp_l%next
 else
  temp_r=>temp_r%previous; temp_l=>temp_l%previous
 end if
 x_idx=temp_l%g_cell%x_idx
 if(x_idx.gt.nxll+3) done=.true.
end do

! call check_list_ac(1,'down  ')

temp_b=>acvv(1)%begin%next%previous
do while(associated(temp_b)) 
 if(temp_b%partner.eq.'single  ') then
  temp_b%n_nxt_dr=>temp_b%next
  temp_b%next%p_nxt_dr=>temp_b
  if(temp_b%next%partner.eq.'next    ') then
   temp_b%next%next%p_nxt_dr=>temp_b
  end if
 else
  if(temp_b%partner.eq.'next    ') then
   temp_b%n_nxt_dr=>temp_b%next%next
  else
   temp_b%n_nxt_dr=>temp_b%next
   temp_b%next%p_nxt_dr=>temp_b
   if(temp_b%next%partner.eq.'next    ') then
    temp_b%next%next%p_nxt_dr=>temp_b
   end if
  end if
 end if
 nullify(temp_b%p_nxt_dr); temp_b=>temp_b%previous
end do

! call check_list_ac(1,'down  ')

temp_b=>acvv(1)%eend%previous%next
do while(associated(temp_b)) 
 if(temp_b%partner.eq.'single  ') then
  temp_b%p_nxt_dr=>temp_b%previous
  temp_b%previous%n_nxt_dr=>temp_b
  if(temp_b%previous%partner.eq.'previous') then
   temp_b%previous%previous%n_nxt_dr=>temp_b
  endif
 else
  if(temp_b%partner.eq.'previous') then
   temp_b%p_nxt_dr=>temp_b%previous%previous
  else
   temp_b%p_nxt_dr=>temp_b%previous
   temp_b%previous%n_nxt_dr=>temp_b
   if(temp_b%previous%partner.eq.'previous') then
    temp_b%previous%previous%n_nxt_dr=>temp_b
   end if      	 	  	   
  end if
 end if	 
 nullify(temp_b%n_nxt_dr); temp_b=>temp_b%next
end do

! call check_list_ac(1,'up    ')

end subroutine material_interface


subroutine insert_boundary_cell(temp,current,position)

implicit none
type(auxi_crit_cell), pointer :: temp, current
character*6, intent(in) :: position

type(adss_info) :: ads

ads=temp%address
nullify(current%previous); nullify(current%next)
acvv(ads%cv_nmb)%total=acvv(ads%cv_nmb)%total+1
current%address%idx=acvv(ads%cv_nmb)%total
select case(position)
 case('before')
! The insert point is in front of the first critical cell in the curve. 
  current%next=>temp
  temp%previous=>current
  acvv(ads%cv_nmb)%begin_boundary%next=>current
 case('after ')
! The insert point is behind the last critical cell in the curve. 
  current%previous=>temp
  temp%next=>current
  acvv(ads%cv_nmb)%eend_boundary%previous=>current
  case default; call error_message
end select

end subroutine insert_boundary_cell


end subroutine Jacobs_PengZ_1


subroutine Jacobs_PengZ_2

implicit none

call smooth_region

call material_interface


contains


subroutine material_interface

implicit none

type(auxi_crit_cell), pointer :: temp_b, temp_e

if(acvv(1)%status.eq.'asleep') return

! Periodic implementation at the head-end.
temp_b=>acvv(1)%begin%next%previous
temp_e=>acvv(1)%eend%previous%previous

do while(associated(temp_b))
 temp_b%p_cell=temp_e%p_cell
 temp_b%g_cell%dis_posi=temp_e%g_cell%dis_posi
 temp_b%g_cell%mdis_posi=temp_e%g_cell%mdis_posi
 temp_b%g_cell%normal=temp_e%g_cell%normal
 temp_b=>temp_b%previous; temp_e=>temp_e%previous
end do
 
! periodic implementation at the tail-end.
temp_b=>acvv(1)%begin%next%next; temp_e=>acvv(1)%eend%previous%next

do while(associated(temp_e))
 temp_e%p_cell=temp_b%p_cell
 temp_e%g_cell%dis_posi=temp_b%g_cell%dis_posi
 temp_e%g_cell%mdis_posi=temp_b%g_cell%mdis_posi
 temp_e%g_cell%normal=temp_b%g_cell%normal
 temp_e=>temp_e%next; temp_b=>temp_b%next
end do

! call check_list_ac(1,'down  ')

end subroutine material_interface


end subroutine Jacobs_PengZ_2


subroutine Jacobs_PengZ_3

implicit none

type(auxi_crit_cell), pointer :: temp_b, temp_e, temp_bn, temp_ep
type(geo_info) :: g_cell_b, g_cell_e
type(phy_info) :: p_cell
real(8) :: xx
real(8), dimension(2) :: nml
integer :: j0, jj

if(acvv(1)%status.eq.'asleep') return

temp_b=>acvv(1)%begin%next
temp_e=>acvv(1)%eend%previous
g_cell_b=temp_b%g_cell; g_cell_e=temp_e%g_cell

if(g_cell_b%g_type.ne.'yy '.or.g_cell_e%g_type.ne.'yy ') call error_message

xx=0.5d0*(temp_b%g_cell%dis_posi(2)+temp_e%g_cell%dis_posi(1))
nml(1)=0.5d0*(temp_b%g_cell%normal(2,1)-temp_e%g_cell%normal(1,1))
nml(2)=0.5d0*(temp_b%g_cell%normal(2,2)+temp_e%g_cell%normal(1,2))
g_cell_b%dis_posi=xx; g_cell_e%dis_posi=xx
g_cell_b%normal(1,1)=-nml(1)
g_cell_b%normal(1,2)=nml(2)
g_cell_b%normal(2,1)=nml(1)
g_cell_b%normal(2,2)=nml(2)
g_cell_e%normal=g_cell_b%normal

p_cell=temp_b%p_cell
p_cell%l_state=0.5d0*(temp_b%p_cell%l_state+temp_e%p_cell%l_state)
p_cell%r_state=0.5d0*(temp_b%p_cell%r_state+temp_e%p_cell%r_state)
p_cell%or_state=0.5d0*(temp_b%p_cell%or_state+temp_e%p_cell%or_state)

temp_b%g_cell=g_cell_b; temp_e%g_cell=g_cell_e
temp_b%p_cell=p_cell; temp_e%p_cell=p_cell

temp_bn=>temp_b%n_nxt_dr
j0=g_cell_b%y_idx; jj=temp_bn%g_cell%y_idx
temp_bn%g_cell%dis_posi(1)=g_cell_b%dis_posi(2)+dfloat(j0-jj)
if(temp_bn%partner.eq.'next    ') then
 jj=temp_bn%next%g_cell%y_idx
 temp_bn%next%g_cell%dis_posi(1)=g_cell_b%dis_posi(2)+dfloat(j0-jj)
end if

temp_ep=>temp_e%p_nxt_dr
j0=g_cell_e%y_idx; jj=temp_ep%g_cell%y_idx
temp_ep%g_cell%dis_posi(2)=g_cell_e%dis_posi(1)+dfloat(j0-jj)
if(temp_ep%partner.eq.'previous') then
 jj=temp_ep%previous%g_cell%y_idx
 temp_ep%previous%g_cell%dis_posi(2)=g_cell_b%dis_posi(1)+dfloat(j0-jj)
end if

end subroutine Jacobs_PengZ_3


subroutine Jacobs_PZ_nonreflecting_1(uum)

implicit none

type(state), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2), intent(inout) :: uum

call smooth_rg_nonreflecting_1(uum)

end subroutine Jacobs_PZ_nonreflecting_1


subroutine Jacobs_PZ_nonreflecting_2(uum,uu1)

implicit none

type(state), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2), intent(in) :: uum, uu1

call smooth_rg_nonreflecting_2(uum,uu1)

end subroutine Jacobs_PZ_nonreflecting_2


subroutine smooth_region

implicit none

integer :: i, j

type(state) :: vvv
real(8), dimension(2) :: normal

nxll_boundary=nxll-3; nx_boundary=nx+3
nyll_boundary=nyll-3; ny_boundary=ny+3

do i=nxll-1,nxll-3,-1
 do j=nyll-3,nyll+ny+2
  if(ggd_cell(i+nx-1,j)%region.eq.'smth  ') then
   uu(i,j)=uu(i+nx-1,j)
   call clean(ggd_cell(i,j))
  end if
 end do
end do

do i=nxll+nx,nxll+nx+2
 do j=nyll-3,nyll+ny+2
  if(ggd_cell(i-nx+1,j)%region.eq.'smth  ') then
   uu(i,j)=uu(i-nx+1,j)
   call clean(ggd_cell(i,j))
  end if
 end do
end do

end subroutine smooth_region


subroutine smooth_rg_nonreflecting_1(uum)

implicit none

type(state), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2), intent(inout) :: uum

real(8), dimension(nxll:nxll+nx-1,0:3) :: rhob, ub, vb, pb, cb 
real(8), dimension(nxll:nxll+nx-1,1:3) :: m1, m2, m3, m4, drhob, dub, dvb, dpb
type(state), dimension(nxll:nxll+nx-1,1:3) :: duub
integer :: i, j

type(auxi_crit_cell), pointer :: temp
type(state), dimension(-3:3) :: ux
type(state), dimension(-2:3) :: data_sector
real(8) :: x_posi, y_posi
type(state) :: fl, fr

! real(8) :: xxx

! Firstly, compute the bottom boundary.

! Produce the primitive variables at the present level.
do i=nxll,nxll+nx-1
 do j=0, 3
  rhob(i,j)=uu(i,nyll-j)%value(1)
  ub(i,j)=x_velocity(uu(i,nyll-j)); vb(i,j)=y_velocity(uu(i,nyll-j))
  pb(i,j)=pressure(uu(i,nyll-j))
  cb(i,j)=sound_speed(uu(i,nyll-j))
 end do
end do

! Compute the m's.
m1=0.0d0; m2=0.0d0; m3=0.0d0; m4=0.0d0
do i=nxll,nxll+nx-1
 do j=1, 3
  if(vb(i,j)-cb(i,j).lt.0.0d0) m1(i,j)=(vb(i,j)-cb(i,j))*(pb(i,j-1)-pb(i,j)-rhob(i,j)*cb(i,j)*(vb(i,j-1)-vb(i,j)))
  if(vb(i,j).lt.0.0d0) then
   m2(i,j)=vb(i,j)*(ub(i,j-1)-ub(i,j)); m3(i,j)=vb(i,j)*(pb(i,j-1)-pb(i,j)-cb(i,j)*cb(i,j)*(rhob(i,j-1)-rhob(i,j)))
  end if
  if(vb(i,j)+cb(i,j).lt.0.0d0) m4(i,j)=(vb(i,j)+cb(i,j))*(pb(i,j-1)-pb(i,j)+rhob(i,j)*cb(i,j)*(vb(i,j-1)-vb(i,j)))
 end do        
end do

! Compute the y-differences of the primitive variables.
do i=nxll,nxll+nx-1
 do j=1, 3
  dpb(i,j)=-0.5d0*(m4(i,j)+m1(i,j))
  dvb(i,j)=-0.5d0*(m4(i,j)-m1(i,j))/rhob(i,j)/cb(i,j)
  dub(i,j)=-m2(i,j)
  drhob(i,j)=(dpb(i,j)+m3(i,j))/cb(i,j)/cb(i,j)
 end do
end do

! Compute the y-differences of the conservative variables.
do i=nxll,nxll+nx-1
 do j=1, 3
  duub(i,j)%value(1)=drhob(i,j)
  duub(i,j)%value(2)=ub(i,j)*drhob(i,j)+rhob(i,j)*dub(i,j)
  duub(i,j)%value(3)=vb(i,j)*drhob(i,j)+rhob(i,j)*dvb(i,j)
  duub(i,j)%value(4)=0.5d0*(ub(i,j)*ub(i,j)+vb(i,j)*vb(i,j))*drhob(i,j)
  duub(i,j)%value(4)=duub(i,j)%value(4)+rhob(i,j)*(ub(i,j)*dub(i,j)+vb(i,j)*dvb(i,j))+dpb(i,j)/(uu(i,j)%gamma-1.0d0)
 end do
end do

! Compute the conservative variables
do i=nxll,nxll+nx-1
 do j=1, 3
  if(ggd_cell(i,nyll-j)%region.ne.'smth  ') call error_message
  nullify(temp)
  call smth_dpr(i,nyll-j,temp,ux,'xx ','      ','ordinary',.true.,3,'forward ')
  x_posi=(dfloat(i)+0.5d0)*h; y_posi=dfloat(nyll-j)*h
  data_sector=ux(-3:2); fl=flux(data_sector,'xx',x_posi,y_posi)
  data_sector=ux(-2:3); fr=flux(data_sector,'xx',x_posi,y_posi)
  uum(i,nyll-j)=uu(i,nyll-j)-r*(fr-fl)+r*duub(i,j)
 end do
end do

! Then, compute the top boundary.

! Produce the primitive variables at the present level.
do i=nxll,nxll+nx-1
 do j=0, 3
  rhob(i,j)=uu(i,nyll+ny-1+j)%value(1)
  ub(i,j)=x_velocity(uu(i,nyll+ny-1+j)); vb(i,j)=y_velocity(uu(i,nyll+ny-1+j))
  pb(i,j)=pressure(uu(i,nyll+ny-1+j))
  cb(i,j)=sound_speed(uu(i,nyll+ny-1+j))
 end do
end do

! Compute the m's.
m1=0.0d0; m2=0.0d0; m3=0.0d0; m4=0.0d0
do i=nxll,nxll+nx-1
 do j=1, 3
  if(vb(i,j)-cb(i,j).gt.0.0d0) m1(i,j)=(vb(i,j)-cb(i,j))*(pb(i,j)-pb(i,j-1)-rhob(i,j)*cb(i,j)*(vb(i,j)-vb(i,j-1)))
  if(vb(i,j).gt.0.0d0) then
   m2(i,j)=vb(i,j)*(ub(i,j)-ub(i,j-1)); m3(i,j)=vb(i,j)*(pb(i,j)-pb(i,j-1)-cb(i,j)*cb(i,j)*(rhob(i,j)-rhob(i,j-1)))
  end if
  if(vb(i,j)+cb(i,j).gt.0.0d0) m4(i,j)=(vb(i,j)+cb(i,j))*(pb(i,j)-pb(i,j-1)+rhob(i,j)*cb(i,j)*(vb(i,j)-vb(i,j-1)))
 end do        
end do

! Compute the y-differences of the primitive variables.
do i=nxll,nxll+nx-1
 do j=1, 3
  dpb(i,j)=-0.5d0*(m4(i,j)+m1(i,j))
  dvb(i,j)=-0.5d0*(m4(i,j)-m1(i,j))/rhob(i,j)/cb(i,j)
  dub(i,j)=-m2(i,j)
  drhob(i,j)=(dpb(i,j)+m3(i,j))/cb(i,j)/cb(i,j)
 end do
end do

! Compute the y-differences of the conservative variables.
do i=nxll,nxll+nx-1
 do j=1, 3
  duub(i,j)%value(1)=drhob(i,j)
  duub(i,j)%value(2)=ub(i,j)*drhob(i,j)+rhob(i,j)*dub(i,j)
  duub(i,j)%value(3)=vb(i,j)*drhob(i,j)+rhob(i,j)*dvb(i,j)
  duub(i,j)%value(4)=0.5d0*(ub(i,j)*ub(i,j)+vb(i,j)*vb(i,j))*drhob(i,j)
  duub(i,j)%value(4)=duub(i,j)%value(4)+rhob(i,j)*(ub(i,j)*dub(i,j)+vb(i,j)*dvb(i,j))+dpb(i,j)/(uu(i,j)%gamma-1.0d0)
 end do
end do

! Compute the conservative variables
do i=nxll,nxll+nx-1
 do j=1, 3
  if(ggd_cell(i,nyll+ny-1+j)%region.ne.'smth  ') call error_message
  nullify(temp)
  call smth_dpr(i,nyll+ny-1+j,temp,ux,'xx ','      ','ordinary',.true.,3,'forward ')
  x_posi=(dfloat(i)+0.5d0)*h; y_posi=dfloat(nyll+ny-1+j)*h
  data_sector=ux(-3:2); fl=flux(data_sector,'xx',x_posi,y_posi)
  data_sector=ux(-2:3); fr=flux(data_sector,'xx',x_posi,y_posi)
  uum(i,nyll+ny-1+j)=uu(i,nyll+ny-1+j)-r*(fr-fl)+r*duub(i,j)
 end do
end do

! xxx=pressure(uum(20,43))
! continue

end subroutine smooth_rg_nonreflecting_1


subroutine smooth_rg_nonreflecting_2(uum,uu1)

implicit none

type(state), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2), intent(in) :: uum, uu1

integer :: i, j

do i=nxll, nxll+nx-1
 do j=1, 3
  uu(i,nyll-j)=0.5d0*(uu1(i,nyll-j)+uum(i,nyll-j))
  uu(i,nyll+ny-1+j)=0.5d0*(uu1(i,nyll+ny-1+j)+uum(i,nyll+ny-1+j))
 end do
end do

end subroutine smooth_rg_nonreflecting_2


end module Jacobs_Peng_Zabusky_boundary