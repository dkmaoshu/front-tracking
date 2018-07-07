module upd_cv_cases_1

use upd_cv_cases
! 'updcvcas.f90'

use node_cells
! 'nd_cls.f90'

implicit none

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: adds_s

public  set_move_2, side_neighbor
private set_move_1  !, update_position, change_neighbor_dp, change_stencil


contains


subroutine set_move_2(t_current)
! Finally determine whether the current critical cell should 
! move, remain or split, and in what way, by taking into account
! neighboring critical cells in the curve list.

implicit none
type(auxi_crit_cell), pointer :: t_current
! The critical cell under concern.

type(auxi_crit_cell), pointer :: t_neighbor, t_nn
! the neighboring critical cell and the neighboring-neighboring 
! critical cell.
logical, dimension(2) :: moved_current, moved_neighbor, moved_nn
type(geo_info) :: g_neighbor, g_nn
character*6 :: happen_current, happen_neighbor, happen_nn

nullify(t_neighbor); nullify(t_nn)
moved_current=.false.; moved_neighbor=.false.; moved_nn=.false.
happen_current='      '; happen_neighbor='      '
happen_nn='      '

! call check_list_ac(2,'down  ')

call set_move_1(t_current,moved_current,happen_current)
call update_position(t_current,1)
call update_position(t_current,2)
if(moved_current(2)) then
 select case(t_current%partner)
  case('single  ','previous') 
   if(associated(t_current%next)) t_neighbor=>t_current%next
  case('next    ')
   if(associated(t_current%next%next)) &
      t_neighbor=>t_current%next%next
  case default; call error_message
 end select
 if(associated(t_neighbor)) then
  g_neighbor=t_neighbor%g_cell
  if(g_neighbor%g_type.eq.'xy ') then
   call set_move_1(t_neighbor,moved_neighbor,happen_neighbor)
   call update_position(t_neighbor,1)
   call update_position(t_neighbor,2)

!    call check_list_ac(3,'down  ')

   if(happen_neighbor.eq.'move  ') then
    if(associated(t_neighbor%next)) then
     t_nn=>t_neighbor%next; g_nn=t_nn%g_cell
	 if(g_nn%g_type.eq.'xy ') then
      call set_move_1(t_nn,moved_nn,happen_nn)
      call update_position(t_nn,1)
      call update_position(t_nn,2)
      call set_move_1(t_neighbor,moved_neighbor,happen_neighbor)
      call update_position(t_neighbor,1)
      call update_position(t_neighbor,2)
     end if
    end if
   end if
   call set_move_1(t_current,moved_current,happen_current)
   call update_position(t_current,1) 
   call update_position(t_current,2)
  end if   
 end if   
end if

! call check_list_ac(2,'down  ')

happen=happen_current


contains


subroutine update_position(temp,number)
! When one of the discontinuity positions in the critical cell was
! changed, the corresponding discontinuity positions in neighboring 
! critical cells should be changed. This subroutine performs this
! operation.

implicit none
type(auxi_crit_cell), pointer :: temp
! The current critical cell, one of the discontinuity positions in
! which has been changed.
integer, intent(in) :: number
! The number of the changed position.

type(auxi_crit_cell), pointer :: t_partner, t_neighbor, tn_partner
! The partner, neighboring, and the neighboring partner critical
! cells.
type(geo_info) :: g_cell, g_partner, g_neighbor, gn_partner
! The corresponding geometrical cells of the above critical cells

! integer, dimension(3) :: stencil_disp

! print*, temp%stencil_disp(3)

g_cell=temp%g_cell
nullify(t_partner); nullify(t_neighbor); nullify(tn_partner)
call clean_up_g_cell(g_partner); call clean_up_g_cell(g_neighbor)
call clean_up_g_cell(gn_partner)

call change_stencil(temp,number)
select case(number)
 case(1)
  select case(temp%partner)
   case('single  ')
  	if(associated(temp%previous)) t_neighbor=>temp%previous
   case('previous')
    t_partner=>temp%previous
    if(associated(t_partner%previous)) t_neighbor=>t_partner%previous
  end select   
  if(associated(t_neighbor)) then
   if(t_neighbor%partner.eq.'previous') tn_partner=>t_neighbor%previous
  end if
 case(2)
  select case(temp%partner)
   case('single  ')
    if(associated(temp%next)) t_neighbor=>temp%next
   case('next    ')
    t_partner=>temp%next
    if(associated(t_partner%next)) t_neighbor=>t_partner%next
  end select   
  if(associated(t_neighbor)) then
   if(t_neighbor%partner.eq.'next    ') tn_partner=>t_neighbor%next
  end if
end select

! call check_list_ac(2,'down  ')
!  print*, t_neighbor%stencil_disp(3)

if(associated(t_neighbor)) then
 g_neighbor=t_neighbor%g_cell
 call change_neighbor_dp(g_cell,g_neighbor,number,'neighbor')
!  stencil_disp=t_neighbor%stencil_disp
 t_neighbor%g_cell=g_neighbor
 call change_stencil(t_neighbor,3-number)
end if
if(associated(t_partner)) then
 g_partner=t_partner%g_cell
 call change_neighbor_dp(g_cell,g_partner,number,'partner ')
 t_partner%g_cell=g_partner
 call change_stencil(t_partner,number)
end if
if(associated(tn_partner)) then
 gn_partner=tn_partner%g_cell
 call change_neighbor_dp(g_cell,gn_partner,number,'neighbor')
 tn_partner%g_cell=gn_partner
 call change_stencil(tn_partner,3-number)
end if

end subroutine update_position


subroutine change_neighbor_dp(g_cell,g_change,number,relation)

implicit none
type(geo_info), intent(in) :: g_cell
type(geo_info), intent(inout) :: g_change
integer, intent(in) :: number
character*8, intent(in) :: relation

integer :: i0, j0, i1, j1
real(8) :: dif

i0=g_cell%x_idx; j0=g_cell%y_idx
i1=g_change%x_idx; j1=g_change%y_idx
select case(g_cell%edge(number))
 case(1,3); dif=dfloat(i0-i1)
 case(2,4); dif=dfloat(j0-j1)
end select
select case(relation)
 case('neighbor')
  g_change%dis_posi(3-number)=g_cell%dis_posi(number)+dif
 case('partner ')
  g_change%dis_posi(number)=g_cell%dis_posi(number)+dif
end select

end subroutine change_neighbor_dp


subroutine change_stencil(tempx,number)
! The change of discontinuity position will cause the stencil in the critical cell to 
! change.

implicit none
type(auxi_crit_cell), pointer :: tempx
integer, intent(in) :: number

type(geo_info) :: gx_cell
real(8), dimension(3,2) :: stencil
integer, dimension(3) :: stencil_disp
real(8), dimension(2) :: pt
integer :: i

! print*

gx_cell=tempx%g_cell; stencil=tempx%stencil; stencil_disp=error_index
do i=1,3
 stencil_disp(i)=tempx%stencil_disp(i)
 
!  print*, stencil_disp(i)

end do

do i=1,3
 if(stencil_disp(i).eq.number) then
  call pick_point(gx_cell,number,pt)
  stencil(i,:)=pt
 end if
end do    

tempx%stencil=stencil

end subroutine change_stencil


end subroutine set_move_2


subroutine side_neighbor(temp,temp_g,head_tail,head_tail_g)
! This subroutine determine on which side of the neighbor critical
! cell the current critical cell is located. It applies only to 
! the case when the current critical cell is either a head or an
! end critical cell of a discontinuity curve and is adjacent to a
! node cell.

type(auxi_crit_cell), pointer :: temp, temp_g
! The pointer pointing to the other critical cell stacked in the
! the same cell with the critical cell pointed by 'temp', the 
! first critical cell must be either a head or an end of a 
! neighbor discontinuity curve pluging in the same node cell.
character*6, intent(in) :: head_tail, head_tail_g
! Indicating whether 'temp' and 'temp_g are a head or an end.

type(curve_plug), pointer :: tempp, tempp_g
type(cv_plug_info) :: plug, plug_g
type(geo_info) :: g_cell
integer :: nd_nmb, cv_nmb
type(adss_info) :: address, address_g

address=temp%address; address_g=temp_g%address

! Find the node cell both the corresponding discontinuity curves
! plug in.
cv_nmb=address%cv_nmb
select case(head_tail)
 case('begin ')
  nd_nmb=acvv(cv_nmb)%begin_end
 case('end   ')
  nd_nmb=acvv(cv_nmb)%end_end
 case default
  print*, 'There must be something wrong.'
  pause
end select

! Find the curve plug corresponding to the discontinuity curve the 
! current critical cell 'temp' belongs to.
g_cell=temp%g_cell
plug%cv_nmb=address%cv_nmb
if(head_tail.eq.'begin ') then
 plug%end_type='begin '; plug%edge=neighbor_edge(g_cell%edge(1))
else
 plug%end_type='end   '; plug%edge=neighbor_edge(g_cell%edge(2))
end if
call visit(nd_nmb,plug,tempp)

! Find the curve plug corresponding to the discontinuity curve the 
! stacked critical cell 'temp_g' belongs to.
g_cell=temp_g%g_cell
plug_g%cv_nmb=address_g%cv_nmb
if(head_tail_g.eq.'begin ') then
 plug_g%end_type='begin '; plug_g%edge=neighbor_edge(g_cell%edge(1))
else
 plug_g%end_type='end   '; plug_g%edge=neighbor_edge(g_cell%edge(2))
end if
call visit(nd_nmb,plug_g,tempp_g)

! Choose one of the pointers 'temp_g%l_stk' and 'temp_g%r_stk'
! pointing to 'temp'.
if(tempp_g%previous%address%idx.eq.tempp%address%idx) then
 select case(plug_g%end_type)
  case('begin '); temp_g%r_stk=address
  case('end   '); temp_g%l_stk=address
 end select
else
 if(tempp_g%next%address%idx.eq.tempp%address%idx) then
  select case(plug_g%end_type)
   case('begin '); temp_g%l_stk=address
   case('end   '); temp_g%r_stk=address
  end select
 else
  print*, 'There must be something wrong.'; pause
 end if
end if

end subroutine side_neighbor


end module upd_cv_cases_1