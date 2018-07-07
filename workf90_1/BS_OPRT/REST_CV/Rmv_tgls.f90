module remove_tangles
! This module removes tangles on the discontinuity curves to be formed and 
! corrects the grid-map accordingly.

use solu_comput
! 'solu_com.f90'

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none

public rm_tg


contains
						

subroutine rm_tg(t_current,g_cell,p_cell,this_side,that_side,if_tangled)
! Remove tangles on a single discontinuity curve and correct 
! the grid-map.

! Tangle is that the current and previous critical cells are 
! located in the same grid cells. To remove it, one needs to merge 
! these neighboring critical cells and modify the related
! quantities accordingly, including the grid-map around the 
! tangles.

implicit none
type(auxi_crit_cell), pointer :: t_current
! The current critical cell.
type(geo_info), intent(in) :: g_cell
! The moved or split geometrical critical cell.
type(phy_info), intent(in) :: p_cell
! The moved or split physical critical cell.
character*6, intent(in) :: this_side,that_side
! The side relative to the tangle.
logical, intent(out) :: if_tangled

type(auxi_crit_cell), pointer :: t_previous, t_xxx
! The previous critical cell.
type(geo_info) :: g_previous, g_xxx
! Geometrical information in the current and next critical cell.
type(phy_info) :: p_previous
! Physical information in the current and next critical cell.
type(state) :: difference, s_value, t_value
integer :: i0, j0, i1, j1

type(state) :: uuu
integer :: iix, jjy, i
logical :: if_facing

if_tangled=.false.
i0=g_cell%x_idx; j0=g_cell%y_idx
if(.not.associated(t_current%previous)) return
t_previous=>t_current%previous
g_previous=t_previous%g_cell; p_previous=t_previous%p_cell
i1=g_previous%x_idx; j1=g_previous%y_idx
! No tangle and return.
if(i0.ne.i1.or.j0.ne.j1) return

if_tangled=.true.
! Remove the tangled critical cell,
!call remove(t_current,'after ','left  ',difference)

! Determine whether the two tangled discontinuity cells are facing to
! each other or not.
do i=1,4
 if(g_cell%point(i).eq.this_side) exit
end do
if(g_previous%point(i).eq.g_cell%point(i)) then
 if_facing=.false.
else
 if_facing=.true.
end if

! Update the physical states.
select case(that_side)
 case('left  ')
  difference=p_cell%or_state-p_cell%l_state
 case('right ')
  difference=p_cell%or_state-p_cell%r_state
end select
p_previous%l_state=0.5d0*(p_cell%l_state+p_previous%l_state)
p_previous%r_state=0.5d0*(p_cell%r_state+p_previous%r_state)
p_previous%or_state=p_previous%or_state+difference

! Update the geometrical information.
g_previous%dis_posi(2)=g_cell%dis_posi(2)
g_previous%edge(2)=g_cell%edge(2)

t_previous%p_cell=p_previous
! Distinguish different cases.
if(g_previous%edge(1).ne.g_cell%edge(2)) then
! The tangle changes the previous critical cell.
 call find_type_from_edges(g_previous)
 call find_points_from_edges(g_previous)
 t_previous%g_cell=g_previous

 t_xxx=>t_current
 if(t_current%partner.eq.'next    ') then
  if(t_current%g_cell%x_idx.ne.i0.or.t_current%g_cell%y_idx.ne.j0) then
   g_xxx=t_current%next%g_cell
   if(g_xxx%x_idx.eq.i0.and.g_xxx%y_idx.eq.j0) then
! Move must be happened and the current one moves to its "next" partner.
    t_xxx=>t_current%next
   end if
  end if
 end if

 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%address.eq.t_xxx%address) then
   ggd_cell(i0,j0)%ccpt(i)%address=t_previous%address
   ggd_cell(i0,j0)%ccpt(i)%side=g_previous%point(i)
  end if 
 end do 

else

!%%%%%%%%%%% patch and correction %%%%%%%%%%%%%% 
 if(if_facing) then
  
!  call remove(t_previous,'before',that_side,difference)
  if_tangled=.false.; return
  
 else
  iix=t_previous%g_cell%x_idx; jjy=t_previous%g_cell%y_idx
  if(this_side.eq.'left  ') then
   uuu=t_previous%p_cell%l_state
  else
   if(this_side.eq.'right ') then
    uuu=t_previous%p_cell%r_state
   else
    call error_message
   end if
  endif
  
  call remove(t_previous,'before',this_side,difference)
  
  uu(iix,jjy)=uuu
!%%%%%%%%%%%% patch and correction %%%%%%%%%%%%%%

 end if  
 
 if(associated(t_previous)) then
  t_previous%p_cell%or_state=t_previous%p_cell%or_state+difference
 end if

end if

end subroutine rm_tg


end module remove_tangles