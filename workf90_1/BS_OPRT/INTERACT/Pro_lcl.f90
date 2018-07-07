module produce_local_smooth
! This module is for production of smooth solution in local in 
! the region between two adjacent discontinuity curves. This 
! smooth solution in local will be used in the handling of 
! interactions of discontinuity curves.

use solution
! 'solution.f90'

use solu_comput
! 'solu_com.f90'

implicit none

type general_crit_cell
 type(critical_cell), pointer :: temp
 type(auxi_crit_cell), pointer :: temp_a
end type general_crit_cell

type(general_crit_cell) :: temp_g, temq_g
! The pointer pointing to the critical cell on the current curve
! and the critical cell on the next curve pluging in the same
! node cell. The critical cells are either regular or auxiliary.
type(critical_cell), pointer :: temp, temq
type(auxi_crit_cell), pointer :: temp_a, temq_a
type(curve_plug), pointer :: temp_p, temq_p
! Pointers pointing to corresponding regular or auxiliary 
! critical cells
type(cv_plug_info) :: plug, plug_1
! Plug-in information corresponding to the two discontinuity
! curves that form the smooth region.
integer :: i0, j0, i1, j1
! The indexes of the node cell and  of the critical cells involved 
! in the production.
type(state), dimension(-1:1,-1:1) :: ulc
! The smooth solution in the local region.
character*6 :: side, side_1
! The side 
real(8),dimension(-3:3,-3:3) :: counter


private i1, j1, first_extension, temp_g, second_extension, &
        produce_smth, counter 
public  produce_smooth_regu, produce_smooth_auxi, i0, j0, &
        side, side_1, uu, plug, plug_1, ulc, temp, temq, &
        temp_a, temq_a, temp_p, temq_p


contains


subroutine produce_smooth_regu
! This subroutine performs the production for regular 
! discontinuity curves.

implicit none

nullify(temp_g%temp_a); temp_g%temp=>temp
nullify(temq_g%temp_a); temq_g%temp=>temq

call produce_smth

end subroutine produce_smooth_regu


subroutine produce_smooth_auxi
! This performs the production for auxiliary discontinuity 
! curves.

implicit none

nullify(temp_g%temp); temp_g%temp_a=>temp_a
nullify(temq_g%temp); temq_g%temp_a=>temq_a

call produce_smth

end subroutine produce_smooth_auxi


subroutine produce_smth
! This subroutine produces the smooth solution in local in the 
! vinicity of critical cell pointed by 'temp' on the assigned
! side.	The local region is bounded by two discontinuity curves
! pluging in the same node cell.

implicit none

type(state), dimension(-3:3,-3:3) :: ux, uy, uuu
integer :: i, j

do i=-1,1; do j=-1,1
 ulc(i,j)=error_data
end do; end do
do i=-3,3; do j=-3,3
 ux(i,j)=error_data; uy(i,j)=error_data; uuu(i,j)=error_data
end do; end do

call first_extension(ux,'xx ')
call first_extension(uy,'yy ')
do i=-3,3; do j=-3,3
 if(ux(i,j).gt.0.9d0*error_data.and. &
    uy(i,j).gt.0.9d0*error_data) then
  uuu(i,j)=0.5d0*(ux(i,j)+uy(i,j))
 else
  if(ux(i,j).gt.0.9d0*error_data) then
   uuu(i,j)=ux(i,j)
  else
   uuu(i,j)=uy(i,j)
  end if
 end if      
end do; end do

call second_extension(uuu)

end subroutine produce_smth


subroutine first_extension(uxy,dir)
! This subrotuine performs the first step of the production of
! smooth solution. It picks the original smooth data and then
! extends it either in x- or y-direction.

implicit none
type(state), dimension(-3:3,-3:3), intent(out) :: uxy
character*3, intent(in) :: dir

! type(geo_info) :: g_cell
! type(phy_info) :: p_cell
type(general_crit_cell) :: temp_nb
integer :: i, j

! logical :: xxx, yyy

! Initialization
do i=-3,3; do j=-3,3
 uxy(i,j)=error_data
end do; end do
counter=0.d0

! First pick smooth data along the current curve, do only with
! non-empty curves.
if(associated(temp_g%temp).or.associated(temp_g%temp_a)) then

! First, pick and extend down curve
 nullify(temp_nb%temp); nullify(temp_nb%temp_a)
 if(associated(temp_g%temp)) temp_nb%temp=>temp_g%temp
 if(associated(temp_g%temp_a)) temp_nb%temp_a=>temp_g%temp_a
! yyy=(associated(temp_g%temp))
! xxx=(associated(temp_nb%temp))

! g_cell=temp_nb%temp_a%g_cell
! p_cell=temp_nb%temp_a%p_cell
 call pick_n_extend(temp_nb,side,'down  ',dir,uxy)

! Then pick and extend up curve.
 if(associated(temp_g%temp)) then
  temp_nb%temp=>temp_g%temp%previous
  nullify(temp_nb%temp_a)
 else
  nullify(temp_nb%temp)
  temp_nb%temp_a=>temp_g%temp_a%previous
 end if
 call pick_n_extend(temp_nb,side,'up    ',dir,uxy)

end if

! Pick and extend along the open egdes of the node cell, the 
! edges on which no discontinuity curve plugs.
call pick_n_extend1(dir,uxy)

! Then pick smooth data along the next curve, do only with 
! non-empty one.
if(associated(temq_g%temp).or.associated(temq_g%temp_a)) then

! First, pick and extend down curve
 nullify(temp_nb%temp); nullify(temp_nb%temp_a)
 if(associated(temq_g%temp)) temp_nb%temp=>temq_g%temp
 if(associated(temq_g%temp_a)) temp_nb%temp_a=>temq_g%temp_a
 call pick_n_extend(temp_nb,side_1,'down  ',dir,uxy)

! Then pick and extend up curve.
 if(associated(temq_g%temp)) then
  temp_nb%temp=>temq_g%temp%previous
 end if
 if(associated(temq_g%temp_a)) then
  temp_nb%temp_a=>temq_g%temp_a%previous
 end if

 call pick_n_extend(temp_nb,side_1,'up    ',dir,uxy)

end if

end subroutine first_extension


subroutine pick_n_extend(temp_nb,side,down_up,dir,uxy)
! This subroutine picks smooth data along the discontinutiy curve
! the critical cell under concern is located in either up- or
! down-curve way.

implicit none
type(general_crit_cell) :: temp_nb
! The pointer pointing to the critical cell to be started with.
character*6, intent(in) :: side
character*6, intent(in) :: down_up
! Indicate whether the pick-n_extend will proceed in up- or down-
! curve way. The possible values for this variable are either
! 'down  ' or 'up    '.
character*3 :: dir
! The direction the pick_n_extend proceeds.
type(state), dimension(-3:3,-3:3), intent(out) :: uxy
! Storage of the produced smooth data.

!type(crit_info) :: crit
type(state), dimension(-3:3) :: ulcc, u_trans
integer :: i=0, j=0
logical :: in_region !, xxx
integer :: kk

! type(geo_info) :: gg_ce

do kk=-3, 3; u_trans(kk)=error_data; end do
do while(associated(temp_nb%temp).or.associated(temp_nb%temp_a))
 ! Determine if the picked critical cell is still in the region.
 in_region=.true.
 if(associated(temp_nb%temp)) then
  i1=temp_nb%temp%g_cell%x_idx; j1=temp_nb%temp%g_cell%y_idx
 else
  i1=temp_nb%temp_a%g_cell%x_idx; j1=temp_nb%temp_a%g_cell%y_idx
 end if

 if(dir.eq.'xx '.and.iabs(j1-j0).gt.3) in_region=.false.
 if(dir.eq.'yy '.and.iabs(i1-i0).gt.3) in_region=.false.

! Pick and extrapolate.
 if(in_region) then
  if(associated(temp_nb%temp)) then
   call smth_dpr(i,j,temp_nb%temp,ulcc,dir,side,'critical',.true.,3,'full    ')
  else
   if(associated(temp_nb%temp_a)) then
    call smth_dpr(i,j,temp_nb%temp_a,ulcc,dir,side,'critical',.true.,3,'full    ')
   else
    print*, 'There must be something wrong.'; pause
   end if
  end if

!   xxx=(associated(temp_nb%temp))

  if(dir.eq.'xx ') then  
   call transmit(ulcc,u_trans,i0,i1)
  else
   call transmit(ulcc,u_trans,j0,j1)
  end if

!   xxx=(associated(temp_nb%temp))

  call fill(u_trans,dir,uxy)
 else
  exit
 end if

! Go to the next critical cell on the curve.
!  gg_ce=temp_nb%temp%g_cell
 select case(down_up)
  case('down  ') 
   if(associated(temp_nb%temp)) then
    temp_nb%temp=>temp_nb%temp%next
   else
    temp_nb%temp_a=>temp_nb%temp_a%next
   end if
  case('up    ')
   if(associated(temp_nb%temp)) then
    temp_nb%temp=>temp_nb%temp%previous
   else
    temp_nb%temp_a=>temp_nb%temp_a%previous
   end if
  case default; print*, 'There must be somthing wrong. '
 end select
end do 

end subroutine pick_n_extend


subroutine pick_n_extend1(dir,uxy)
! This subroutine picks smooth data along open edges of the node
! cell, the edges with no discontinuity curve pluging in, between
! the two neighbor discontinuity curves.

implicit none
character*3 :: dir
! The direction the pick_n_extend proceeds.
type(state), dimension(-3:3,-3:3), intent(out) :: uxy
! Storage of the produced smooth data.

type(critical_cell), pointer :: temp
type(auxi_crit_cell), pointer :: temp_a
type(corner_pointer) :: ccpt
type(state), dimension(-3:3) :: ulcc, u_trans
character*6 :: side
integer :: edge_o, kk

do kk=-3,3; u_trans(kk)=error_data; end do

if(plug%edge.eq.plug_1%edge) return
edge_o=cycle_4(plug%edge+1)
do while(edge_o.lt.plug_1%edge)

! Select the corresponding grid cell. 
 select case(edge_o)
  case(1); i1=i0; j1=j0-1
  case(2); i1=i0+1; j1=j0
  case(3); i1=i0; j1=j0+1
  case(4); i1=i1-1; j1=j0
 end select

 nullify(temp); nullify(temp_a)
 side='      '
! Pick and extend.
 ccpt=ggd_cell(i1,j1)%ccpt(edge_o)
 if(ccpt%address%cv_nmb.eq.0) then
  if(associated(temp_g%temp)) then
   call smth_dpr(i1,j1,temp,ulcc,dir,side,'ordinary',.true.,3,'full    ')
  else
   call smth_dpr(i1,j1,temp_a,ulcc,dir,side,'ordinary',.true.,3,'full    ')
  end if
 else
  if(associated(temp_g%temp)) then
   call visit(ccpt%address,temp)
   side=ccpt%side
   call smth_dpr(i1,j1,temp,ulcc,dir,side,'critical',.true.,3,'full    ')
  else
   call visit(ccpt%address,temp_a)
   side=ccpt%side
   call smth_dpr(i1,j1,temp_a,ulcc,dir,side,'critical',.true.,3,'full    ')
  endif
 end if

 if(dir.eq.'xx ') then  
  call transmit(ulcc,u_trans,i0,i1)
 else
  call transmit(ulcc,u_trans,j0,j1)
 end if
 call fill(u_trans,dir,uxy)

 edge_o=edge_o+1
end do 

end subroutine pick_n_extend1


subroutine transmit(uv,u_trans,idx0,idx1)
! This subroutine transmit the smooth data obtained by 'smth_dpr'
! to the correct row or column in the local region.

implicit none
type(state), dimension(-3:3), intent(in) :: uv
type(state), dimension(-3:3), intent(out) :: u_trans
integer, intent(in) :: idx0, idx1

integer :: i

do i=-3,3; u_trans(i)=error_data; end do
do i=max0(0,idx1-idx0)-3,min0(0,idx1-idx0)+3 
u_trans(i)=uv(i+idx0-idx1)
end do

end subroutine transmit


subroutine fill(u_trans,dir,uxy)
! This subroutine fills the smooth data into the local region. 

implicit none
type(state), dimension(-3:3), intent(in) :: u_trans
!integer, dimension(-3:3,-3:3), intent(inout) :: counter
character*3, intent(in) :: dir
type(state), dimension(-3:3,-3:3) :: uxy

integer :: i, j

if(dir.eq.'xx ') then
 do i=-3,3
  if(u_trans(i).gt.0.9d0*error_data) then
   uxy(i,j1-j0)=counter(i,j1-j0)/(counter(i,j1-j0)+1.0d0)* &
                uxy(i,j1-j0)
   uxy(i,j1-j0)=uxy(i,j1-j0)+1.d0/(counter(i,j1-j0)+1.0d0)*u_trans(i)
   counter(i,j1-j0)=counter(i,j1-j0)+1.d0
  end if
 end do
else
 do j=-3,3
  if(u_trans(j).gt.0.9d0*error_data) then
   uxy(i1-i0,j)=counter(i1-i0,j)/(counter(i1-i0,j)+1.0d0)* &
                uxy(i1-i0,j)
   uxy(i1-i0,j)=uxy(i1-i0,j)+1.d0/(counter(i1-i0,j)+1.0d0)*u_trans(j)
   counter(i1-i0,j)=counter(i1-i0,j)+1.d0
  end if
 end do
end if

end subroutine fill


subroutine second_extension(uuu)
! This subroutine performs the second step of the production of 
! smooth solution. Wherever the data of smooth solution has not
! been produced, the subroutine fills it with extrapolation 
! data.

type(state), dimension(-3:3,-3:3), intent(in) :: uuu

type(state), dimension(-1:1,-1:1) :: ulcc
type(state), dimension(-3:3) :: uxy
type(state) :: vx, vy
integer :: i, j, ii, jj

do i=-1,1; do j=-1,1
 ulcc(i,j)=uuu(i,j); ulc(i,j)=uuu(i,j)
end do; end do

do i=-1,1; do j=-1,1
 if(ulcc(i,j).lt.0.9d0*error_data) then
  vx=error_data; vy=error_data

! Compute extrapolation datum on x-direction. 
  do ii=-3,3; uxy(ii)=error_data; end do 
  do ii=-3,3
   if(uuu(ii,j).gt.0.9d0*error_data) then
    uxy(ii)=uuu(ii,j)
   end if
  end do
  if(uxy(0).gt.0.9d0*error_data) then
   call data_fill(uxy,3,'full    ')
   vx=uxy(i)
  end if

! Compute extrapolation datum on y-direction. 
  do jj=-3,3; uxy(jj)=error_data; end do
  do jj=-3,3
   if(uuu(i,jj).gt.0.9d0*error_data) then
    uxy(jj)=uuu(i,jj)
   end if
  end do
  if(uxy(0).gt.0.9d0*error_data) then
   call data_fill(uxy,3,'full    ')
   vy=uxy(j)
  end if

! Average the extrapolation data on two directions.
  if(vx.gt.0.9d0*error_data.and.vy.gt.0.9d0*error_data) then
   ulc(i,j)=0.5d0*(vx+vy)
  else
   if(vx.gt.0.9d0*error_data) then
    ulc(i,j)=vx
   else
    if(vy.gt.0.9d0*error_data) then
     ulc(i,j)=vy
!    else
!     print*, 'There must be something wrong.'
    end if
   end if
  end if
 
 end if 
end do; end do

end subroutine second_extension

end module produce_local_smooth