module Haas_Sturtevant_1d_boundary

use solu_comput
! 'solu_comp.f90'

implicit none

public  Haas_Sturtevant_1d_1, Haas_Sturtevant_1d_2, Haas_Sturtevant_1d_3
private smooth_region
 

contains


subroutine Haas_Sturtevant_1d_1

implicit none

call smooth_region

call material_interface


contains


subroutine material_interface

implicit none

type(auxi_crit_cell), pointer :: temp_b, temp_t, current, temp_bd
integer :: i, j, x_idx, y_idx, x_ii, y_jj, ii, jj
logical :: done, head_bottom

do ii=1, 2
! Determine whether the head of the interface is on the top or bottom.
 i=acvv(ii)%begin%next%g_cell%y_idx
 if(i.eq.nyll) then
  head_bottom=.true.
 else
  if(i.eq.nyll+ny-1) then
   head_bottom=.false.
  else
   call error_message
  end if
 end if
   
! In-out extension on the bottom.
 if(head_bottom) then
  temp_b=>acvv(ii)%begin%next
 else
  temp_b=>acvv(ii)%eend%previous
 end if
  
 do jj=1, 3
  allocate(current)
  current=temp_b
  current%g_cell%y_idx=nyll-jj
  if(head_bottom) then
   call insert_boundary_cell(temp_b,current,'before')
  else
   call insert_boundary_cell(temp_b,current,'after ')
  end if
  x_idx=current%g_cell%x_idx; y_idx=current%g_cell%y_idx
  ggd_cell(x_idx,y_idx)%region='crit  '
  x_ii=temp_b%g_cell%x_idx; y_jj=temp_b%g_cell%y_idx
  do j=1, 4
   ggd_cell(x_idx,y_idx)%ccpt(j)%address=current%address
   ggd_cell(x_idx,y_idx)%ccpt(j)%side=ggd_cell(x_ii,y_jj)%ccpt(j)%side
  end do
  if(head_bottom) then
   temp_b=>temp_b%previous
  else
   temp_b=>temp_b%next
  end if
 end do

! In-out extension on the top.
 if(head_bottom) then
  temp_t=>acvv(ii)%eend%previous
 else
  temp_t=>acvv(ii)%begin%next
 end if

 do jj=1, 3
  allocate(current)
  current=temp_t
  current%g_cell%y_idx=nyll+ny+jj-1
  if(head_bottom) then
   call insert_boundary_cell(temp_t,current,'after ')
  else
   call insert_boundary_cell(temp_t,current,'before')
  end if
  x_idx=current%g_cell%x_idx; y_idx=current%g_cell%y_idx
  ggd_cell(x_idx,y_idx)%region='crit  '
  x_ii=temp_t%g_cell%x_idx; y_jj=temp_t%g_cell%y_idx
  do j=1, 4
   ggd_cell(x_idx,y_idx)%ccpt(j)%address=current%address
   ggd_cell(x_idx,y_idx)%ccpt(j)%side=ggd_cell(x_ii,y_jj)%ccpt(j)%side
  end do
  if(head_bottom) then
   temp_t=>temp_t%next
  else
   temp_t=>temp_t%previous
  end if
 end do

! call check_list_ac(1,'down  ')

 temp_bd=>acvv(ii)%begin%next%previous
 do while(associated(temp_bd)) 
  if(temp_bd%partner.eq.'single  ') then
   temp_bd%n_nxt_dr=>temp_bd%next
   temp_bd%next%p_nxt_dr=>temp_bd
   if(temp_bd%next%partner.eq.'next    ') then
    temp_bd%next%next%p_nxt_dr=>temp_bd
   end if
  else
   if(temp_bd%partner.eq.'next    ') then
    temp_bd%n_nxt_dr=>temp_bd%next%next
   else
    temp_bd%n_nxt_dr=>temp_bd%next
    temp_bd%next%p_nxt_dr=>temp_bd
    if(temp_bd%next%partner.eq.'next    ') then
     temp_bd%next%next%p_nxt_dr=>temp_bd
    end if
   end if
  end if
  nullify(temp_b%p_nxt_dr); temp_bd=>temp_bd%previous
 end do

! call check_list_ac(1,'down  ')

 temp_bd=>acvv(ii)%eend%previous%next
 do while(associated(temp_bd)) 
  if(temp_bd%partner.eq.'single  ') then
   temp_bd%p_nxt_dr=>temp_bd%previous
   temp_bd%previous%n_nxt_dr=>temp_bd
   if(temp_bd%previous%partner.eq.'previous') then
    temp_bd%previous%previous%n_nxt_dr=>temp_bd
   endif
  else
   if(temp_bd%partner.eq.'previous') then
    temp_bd%p_nxt_dr=>temp_bd%previous%previous
   else
    temp_bd%p_nxt_dr=>temp_bd%previous
    temp_bd%previous%n_nxt_dr=>temp_bd
    if(temp_bd%previous%partner.eq.'previous') then
     temp_bd%previous%previous%n_nxt_dr=>temp_bd
    end if      	 	  	   
   end if
  end if	 
  nullify(temp_b%n_nxt_dr); temp_bd=>temp_bd%next
 end do

!  call check_list_ac(1,'down  ')

end do

! call check_list_ac(2,'down  ')

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


end subroutine Haas_Sturtevant_1d_1


subroutine Haas_Sturtevant_1d_2

implicit none

call smooth_region

call material_interface


contains


subroutine material_interface

implicit none

type(auxi_crit_cell), pointer :: temp_b, temp_e, temp_bb, temp_ee

integer :: ii


! Inout implementation at the head-end.
do ii=1, 2

 temp_bb=>acvv(ii)%begin%next
 temp_bb%g_cell%dis_posi(1)=temp_bb%g_cell%dis_posi(2)
 temp_b=>acvv(ii)%begin%next%previous

 do while(associated(temp_b))
  temp_b%p_cell=temp_bb%p_cell
  temp_b%g_cell%dis_posi=temp_bb%g_cell%dis_posi
  temp_b%g_cell%normal=temp_bb%g_cell%normal
  temp_b=>temp_b%previous
 end do
 
! Inout implementation at the tail-end.
 temp_ee=>acvv(ii)%eend%previous
 temp_ee%g_cell%dis_posi(2)=temp_ee%g_cell%dis_posi(1)
 temp_e=>acvv(ii)%eend%previous%next

 do while(associated(temp_e))
  temp_e%p_cell=temp_ee%p_cell
  temp_e%g_cell%dis_posi=temp_ee%g_cell%dis_posi
  temp_e%g_cell%normal=temp_ee%g_cell%normal
  temp_e=>temp_e%next
 end do

end do

! call check_list_ac(1,'down  ')

end subroutine material_interface


end subroutine Haas_sturtevant_1d_2


subroutine Haas_Sturtevant_1d_3

implicit none

!type(auxi_crit_cell), pointer :: temp_b, temp_e, temp_bn, temp_ep
!type(geo_info) :: g_cell_b, g_cell_e
!type(phy_info) :: p_cell
!real(8) :: xx
!real(8), dimension(2) :: nml
!integer :: j0, jj

!if(acvv(1)%status.eq.'asleep') return

!temp_b=>acvv(1)%begin%next
!temp_e=>acvv(1)%eend%previous
!g_cell_b=temp_b%g_cell; g_cell_e=temp_e%g_cell

!if(g_cell_b%g_type.ne.'yy '.or.g_cell_e%g_type.ne.'yy ') call error_message

!xx=0.5d0*(temp_b%g_cell%dis_posi(2)+temp_e%g_cell%dis_posi(1))
!nml(1)=0.5d0*(temp_b%g_cell%normal(2,1)-temp_e%g_cell%normal(1,1))
!nml(2)=0.5d0*(temp_b%g_cell%normal(2,2)+temp_e%g_cell%normal(1,2))
!g_cell_b%dis_posi=xx; g_cell_e%dis_posi=xx
!g_cell_b%normal(1,1)=-nml(1)
!g_cell_b%normal(1,2)=nml(2)
!g_cell_b%normal(2,1)=nml(1)
!g_cell_b%normal(2,2)=nml(2)
!g_cell_e%normal=g_cell_b%normal

!p_cell=temp_b%p_cell
!p_cell%l_state=0.5d0*(temp_b%p_cell%l_state+temp_e%p_cell%l_state)
!p_cell%r_state=0.5d0*(temp_b%p_cell%r_state+temp_e%p_cell%r_state)
!p_cell%or_state=0.5d0*(temp_b%p_cell%or_state+temp_e%p_cell%or_state)

!temp_b%g_cell=g_cell_b; temp_e%g_cell=g_cell_e
!temp_b%p_cell=p_cell; temp_e%p_cell=p_cell

!temp_bn=>temp_b%n_nxt_dr
!j0=g_cell_b%y_idx; jj=temp_bn%g_cell%y_idx
!temp_bn%g_cell%dis_posi(1)=g_cell_b%dis_posi(2)+dfloat(j0-jj)
!if(temp_bn%partner.eq.'next    ') then
! jj=temp_bn%next%g_cell%y_idx
! temp_bn%next%g_cell%dis_posi(1)=g_cell_b%dis_posi(2)+dfloat(j0-jj)
!end if

!temp_ep=>temp_e%p_nxt_dr
!j0=g_cell_e%y_idx; jj=temp_ep%g_cell%y_idx
!temp_ep%g_cell%dis_posi(2)=g_cell_e%dis_posi(1)+dfloat(j0-jj)
!if(temp_ep%partner.eq.'previous') then
! jj=temp_ep%previous%g_cell%y_idx
! temp_ep%previous%g_cell%dis_posi(2)=g_cell_b%dis_posi(1)+dfloat(j0-jj)
!end if

end subroutine Haas_Sturtevant_1d_3


subroutine smooth_region

implicit none

integer :: i, j

nxll_boundary=nxll-3; nx_boundary=nx+3
nyll_boundary=nyll-3; ny_boundary=ny+3

! The left is an "in-out" boundary.
do i=nxll-1,nxll-3,-1
 do j=nyll,nyll+ny-1
  uu(i,j)=uu(nxll,j)
  call clean(ggd_cell(i,j))
 end do
end do

! The right is an "in-out" boundary.
do i=nxll+nx,nxll+nx+2
 do j=nyll,nyll+ny-1
  uu(i,j)=uu(nxll+nx-1,j)
  call clean(ggd_cell(i,j))
 end do
end do

! The bottom is an "in-out" boundary.
do i=nxll, nxll+nx-1
 do j=nyll-3, nyll-1
  if(ggd_cell(i,nyll)%region.eq.'smth  ') then
   uu(i,j)=uu(i,nyll)
   call clean(ggd_cell(i,j))
  end if
 end do
end do

! The top is an "in-out" boundary.
do i=nxll, nxll+nx-1
 do j=nyll+ny,nyll+ny+2
  if(ggd_cell(i,nyll+ny-1)%region.eq.'smth  ') then
   uu(i,j)=uu(i,nyll+ny-1)
   call clean(ggd_cell(i,j))
  end if
 end do
end do

end subroutine smooth_region


end module Haas_Sturtevant_1d_boundary