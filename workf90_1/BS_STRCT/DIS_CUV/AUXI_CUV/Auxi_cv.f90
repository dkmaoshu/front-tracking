module auxiliary_discontinuity_curves

use grid, cvn=>curves_number
! 'grid.f90'

use auxiliary_critical_cells
! 'auxi_cv.f90'

implicit none

type auxi_discv
 type(auxi_crit_cell), pointer :: begin, eend
! The doubly-linked list of critical cells of a discontinuity 
! curve 'cv' has the following structure: The begin node 
! 'cv%begin' points to the first critical cell 'cv%begin%next'
! of the curve and the end node 'cv%eend' points to the last
! criticel cell 'cv%eend%previous'; however, both the previous
! of the first critical cell and the next of the last one,
! 'cv%begin%next%previous' and 'cv%eend%previous%next' are
! null.
 type(auxi_crit_cell), pointer :: begin_boundary, eend_boundary
! Boundary conditions may be imposed on the tracked discontinuity
! curves, and thus the curves need to be extended out of the solution
! region. These two curve heads are used to point the extension part
! of the curve
 character*6 :: status
! The possible values for variable 'status' are 'asleep', 
! 'awake' and 'yawn  '.
 character*8 :: cv_type
! The possible values for variable 'cv_type' are 'dblink' and 
! 'circular'.
 integer :: begin_end, end_end
! The numbers of the nodes the two ends end; if zero, the end
! ends in infinitive boundary.
 integer :: total
! The total number of the critical cells on the curve.
 integer :: wave
end type auxi_discv

type(auxi_discv), dimension(cvn) :: acv
type(auxi_discv), dimension(:), allocatable :: acvv

interface visit
 module procedure visit_ac
end interface

interface deletee
 module procedure delete_ac, delete_ac_node
end interface

interface insert
 module procedure insert_ac
end interface

interface assignment(=)
 module procedure assign_curve
end interface

private acv, visit_ac, delete_ac, delete_ac_node, insert_ac, &
        assign_curve
public  acvv, give_acv, get_acv, check_list_ac, visit, deletee, &
        insert, creat_acv, hidden_state_check_a, side_area_length_partnered


contains


subroutine give_acv

implicit none
integer :: i

do i=1,cvn
 acvv(i)%status=acv(i)%status
 acvv(i)%cv_type=acv(i)%cv_type
 acvv(i)%begin_end=acv(i)%begin_end; acvv(i)%end_end=acv(i)%end_end
 acvv(i)%total=acv(i)%total; acvv(i)%wave=acv(i)%wave
 nullify(acvv(i)%begin); nullify(acvv(i)%eend)
 if(acv(i)%status.eq.'awake '.or.acv(i)%status.eq.'yawn  ') then
  acvv(i)%begin=>acv(i)%begin; acvv(i)%eend=>acv(i)%eend
  acvv(i)%begin_boundary=>acv(i)%begin_boundary 
  acvv(i)%eend_boundary=>acv(i)%eend_boundary
 end if
end do

end subroutine give_acv


subroutine get_acv

implicit none 
integer :: i

do i=1,cvn
 acv(i)%status=acvv(i)%status
 acv(i)%cv_type=acvv(i)%cv_type
 acv(i)%begin_end=acvv(i)%begin_end; acv(i)%end_end=acvv(i)%end_end
 acv(i)%total=acvv(i)%total; acv(i)%wave=acvv(i)%wave
 nullify(acv(i)%begin); nullify(acv(i)%eend)
 if(acvv(i)%status.eq.'awake '.or.acvv(i)%status.eq.'yawn  ') then
  acv(i)%begin=>acvv(i)%begin; acv(i)%eend=>acvv(i)%eend
  acv(i)%begin_boundary=>acvv(i)%begin_boundary 
  acv(i)%eend_boundary=>acvv(i)%eend_boundary
 end if
end do

end subroutine get_acv


subroutine check_list_ac(cv_nmb,up_down)

implicit none

character*6, intent(in) :: up_down
integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: temp, temp_g
type(geo_info) :: gc_c, gc_p, gc_n, gc_l, gc_r, gc_pnd, gc_nnd
type(phy_info) :: pc_c, pc_p, pc_n, pc_l, pc_r, pc_pnd, pc_nnd
type(phy_info) :: phy_state
type(comp_info) :: c_cell
type(flux_info) :: f_cell
type(adss_info) :: address
integer :: head_mark
logical :: head_switch
character*8 :: partner
real(8), dimension(2) :: normal, pt1, pt2

! type(auxi_crit_cell), pointer :: ttxx

! ttxx=>acvv(1)%begin%next%previous
! ttxx=>acvv(1)%eend%previous%next

if(acvv(cv_nmb)%status.eq.'asleep') then
 print*, 'The curve is still not waked up yet.'; return
end if

if(up_down.eq.'up') then
 if(associated(acvv(cv_nmb)%eend%previous)) then
  temp=>acvv(cv_nmb)%eend_boundary%previous
  head_mark=temp%address%idx; head_switch=.true.
 else
  print*, 'The auxiliary curve list is empty.'; return
 end if
else
 if(associated(acvv(cv_nmb)%begin%next)) then
  temp=>acvv(cv_nmb)%begin_boundary%next
  head_mark=temp%address%idx; head_switch=.true.
 else
  print*, 'The auxiliary curve list is empty.'; return
 end if
end if  

do while(associated(temp).and.head_switch)
 call clean_up_g_cell(gc_c); call clean_up_p_cell(pc_c)
 call clean_up_g_cell(gc_n); call clean_up_p_cell(pc_n)
 call clean_up_g_cell(gc_p); call clean_up_p_cell(pc_p)
 call clean_up_g_cell(gc_l); call clean_up_p_cell(pc_l)
 call clean_up_g_cell(gc_r); call clean_up_p_cell(pc_r)
 call clean_up_g_cell(gc_pnd); call clean_up_p_cell(pc_pnd)
 call clean_up_g_cell(gc_nnd); call clean_up_p_cell(pc_nnd)
 gc_c=temp%g_cell; pc_c=temp%p_cell; address=temp%address
 c_cell=temp%c_cell; f_cell=temp%f_cell; partner=temp%partner
 if(associated(temp%previous)) then
  gc_p=temp%previous%g_cell; pc_p=temp%previous%p_cell
 end if
 if(associated(temp%next)) then
  gc_n=temp%next%g_cell; pc_n=temp%next%p_cell
 end if
 if(temp%l_stk%cv_nmb.gt.0) then
  call visit(temp%l_stk,temp_g)
  gc_l=temp_g%g_cell; pc_l=temp_g%p_cell
 endif
 if(temp%r_stk%cv_nmb.gt.0) then
  call visit(temp%r_stk,temp_g)
  gc_r=temp_g%g_cell; pc_r=temp_g%p_cell
 end if
 if(associated(temp%p_nxt_dr)) then
  gc_pnd=temp%p_nxt_dr%g_cell; pc_pnd=temp%p_nxt_dr%p_cell
 end if
 if(associated(temp%n_nxt_dr)) then
  gc_nnd=temp%n_nxt_dr%g_cell; pc_nnd=temp%n_nxt_dr%p_cell
 end if

!  call pick_point(gc_c,1,pt1); call pick_point(gc_c,2,pt2)
!  normal=normal_of_line(pt1,pt2)

  normal=(/1.0d0,0.0d0/)

 call find_physical_state(pc_c%l_state,normal,phy_state%l_state)
 call find_physical_state(pc_c%r_state,normal,phy_state%r_state)

 if(pc_c%or_state%gamma.gt.0.0d0) call error_message

 if(up_down.eq.'up') then
  temp=>temp%previous
 else
  temp=>temp%next
 end if
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if

end do

end subroutine check_list_ac


subroutine check_partner_conserv_ac(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: temp, temp_g
type(geo_info) :: gc_c, gc_p, gc_n, gc_l, gc_r
type(phy_info) :: pc_c, pc_p, pc_n, pc_l, pc_r
type(adss_info) :: address
integer :: head_mark
logical :: head_switch
character*8 :: partner
real(8), dimension(4) :: xxx, yyy
integer :: ixy

if(acvv(cv_nmb)%status.eq.'asleep') then
 print*, 'The curve is still not waked up yet.'; return
end if

if(associated(acvv(cv_nmb)%begin%next)) then
 temp=>acvv(cv_nmb)%begin_boundary%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The auxiliary curve list is empty.'; return
end if  

do while(associated(temp).and.head_switch)
 call clean_up_g_cell(gc_c); call clean_up_p_cell(pc_c)
 call clean_up_g_cell(gc_n); call clean_up_p_cell(pc_n)
 call clean_up_g_cell(gc_p); call clean_up_p_cell(pc_p)
 call clean_up_g_cell(gc_l); call clean_up_p_cell(pc_l)
 call clean_up_g_cell(gc_r); call clean_up_p_cell(pc_r)
 
 gc_c=temp%g_cell; pc_c=temp%p_cell; address=temp%address
 if(associated(temp%previous)) then
  gc_p=temp%previous%g_cell; pc_p=temp%previous%p_cell
 end if
 if(associated(temp%next)) then
  gc_n=temp%next%g_cell; pc_n=temp%next%p_cell
 end if
 if(temp%l_stk%cv_nmb.gt.0) then
  call visit(temp%l_stk,temp_g)
  gc_l=temp_g%g_cell; pc_l=temp_g%p_cell
 endif
 if(temp%r_stk%cv_nmb.gt.0) then
  call visit(temp%r_stk,temp_g)
  gc_r=temp_g%g_cell; pc_r=temp_g%p_cell
 end if
   
 if(temp%partner.ne.'single   ') then
   
  if(temp%partner.eq.'previous') then
   if(gc_p%x_idx-gc_c%x_idx+gc_p%y_idx-gc_c%y_idx.gt.0) then
    if(gc_c%point(3).eq.'left  ') then
     xxx=pc_c%or_state%value+pc_p%l_state%value
     yyy=pc_p%or_state%value+pc_c%r_state%value
    else
     xxx=pc_c%or_state%value+pc_p%r_state%value
     yyy=pc_p%or_state%value+pc_c%l_state%value
    end if
   else
    if(gc_c%point(1).eq.'left  ') then
     xxx=pc_c%or_state%value+pc_p%l_state%value
     yyy=pc_p%or_state%value+pc_c%r_state%value
    else
     xxx=pc_c%or_state%value+pc_p%r_state%value
     yyy=pc_p%or_state%value+pc_c%l_state%value
    end if
   end if
  else
   if(gc_n%x_idx-gc_c%x_idx+gc_n%y_idx-gc_c%y_idx.gt.0) then
    if(gc_c%point(3).eq.'left  ') then
     xxx=pc_c%or_state%value+pc_n%l_state%value
     yyy=pc_n%or_state%value+pc_c%r_state%value
    else
     xxx=pc_c%or_state%value+pc_n%r_state%value
     yyy=pc_n%or_state%value+pc_c%l_state%value
    end if
   else
    if(gc_c%point(1).eq.'left  ') then
     xxx=pc_c%or_state%value+pc_n%l_state%value
     yyy=pc_n%or_state%value+pc_c%r_state%value
    else
     xxx=pc_c%or_state%value+pc_n%r_state%value
     yyy=pc_n%or_state%value+pc_c%l_state%value
    end if
   end if   	 
  end if

  if(temp%l_stk%cv_nmb.eq.0.and.temp%r_stk%cv_nmb.eq.0) then 
   if(temp%partner.eq.'previous') then
    if(temp%previous%l_stk%cv_nmb.eq.0.and.temp%previous%r_stk%cv_nmb.eq.0) then
     do ixy=1, 4
      if(dabs(xxx(ixy)-yyy(ixy)).gt.1.0d-12)  call error_message
     end do
    end if
   else
    if(temp%next%l_stk%cv_nmb.eq.0.and.temp%next%r_stk%cv_nmb.eq.0) then
     do ixy=1, 4
      if(dabs(xxx(ixy)-yyy(ixy)).gt.1.0d-12)  call error_message
     end do
    end if
   end if
  end if

 end if 

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if

end do

end subroutine check_partner_conserv_ac


subroutine visit_ac(adss,temp_g)
! The subroutine for visit of auxiliary critical cells by their 
! addresses. Variable 'adss' is the address of the auxiliary 
! critical cell to be visited, and variable 'temp_g' is the pointer 
! to the quest auxiliary critical cell containing information 
! passed in or out.

implicit none

type(adss_info), intent(in)	:: adss
type(auxi_crit_cell), pointer :: temp_g

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(flux_info) :: f_cell
type(comp_info) :: c_cell
character*8 :: partner
type(adss_info) :: l_stk, r_stk, address
integer :: head_mark
logical :: head_switch

if(acvv(adss%cv_nmb)%status.eq.'asleep ') then
 print*, 'The curve list does not exist.'; return
end if

if(associated(acvv(adss%cv_nmb)%begin%next)) then
 temp=>acvv(adss%cv_nmb)%begin_boundary%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if

do while(associated(temp).and.head_switch) 
 address=temp%address
 if(address%idx.eq.adss%idx) then
  temp_g=>temp
  g_cell=temp_g%g_cell; p_cell=temp_g%p_cell
  f_cell=temp_g%f_cell; c_cell=temp_g%c_cell
  partner=temp_g%partner
  l_stk=temp_g%l_stk; r_stk=temp_g%r_stk
  exit
 endif       
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  print*, 'No match in the visit.'
  call error_message
  exit
 endif

end do   

end subroutine visit_ac


subroutine hidden_state_check_a

implicit none

integer :: cv_nmb
type(auxi_crit_cell), pointer :: temp, temp_n
type(geo_info) :: g_cell, g_cell_n
type(phy_info) :: p_cell, p_cell_n 
integer :: head_mark
logical :: head_switch
real(8) :: difference

do cv_nmb=1,cvn
 if(acvv(cv_nmb)%status.eq.'asleep') return
  
 nullify(temp); nullify(temp_n)
  
 if(associated(acvv(cv_nmb)%begin%next)) then
  temp=>acvv(cv_nmb)%begin%next
  head_mark=temp%address%idx; head_switch=.true.
 else
  return
 end if  
    
 do while(associated(temp).and.head_switch)
  
  g_cell=temp%g_cell; p_cell=temp%p_cell
  
  if(temp%l_stk%cv_nmb.ne.0) then
   call visit(temp%l_stk,temp_n)
   g_cell_n=temp_n%g_cell; p_cell_n=temp_n%p_cell
   if(temp_n%l_stk.eq.temp%address) then
    difference=dabs(p_cell%l_state-p_cell_n%l_state)
   else
    difference=dabs(p_cell%l_state-p_cell_n%r_state)
   end if
   if(difference.gt.1.0d-9) call error_message
  end if
  if(temp%r_stk%cv_nmb.ne.0) then
   call visit(temp%r_stk,temp_n)
   g_cell_n=temp_n%g_cell; p_cell_n=temp_n%p_cell
   if(temp_n%l_stk.eq.temp%address) then
    difference=dabs(p_cell%r_state-p_cell_n%l_state)
   else
    difference=dabs(p_cell%r_state-p_cell_n%r_state)
   end if
   if(difference.gt.1.0d-9) call error_message
  end if                  
   
  temp=>temp%next
  if(associated(temp)) then
   head_switch=(temp%address%idx.ne.head_mark)
  else
   exit
  end if
  
 end do
  
end do

end subroutine hidden_state_check_a


subroutine delete_ac(acv)
! Delete the whole curve.

implicit none
type(auxi_discv), intent(inout) :: acv

type(auxi_crit_cell), pointer :: temp, current

! type(geo_info) :: g_cell
! logical :: xxx 

if(acv%status.eq.'awake '.or.acv%status.eq.'yawn  ') then
 acv%status='asleep'
 acv%cv_type=' '
 acv%begin_end=-1000; acv%end_end=-1000
 acv%total=-1000

 if(associated(acv%begin%next)) then
  acv%begin_boundary%next%previous=>acv%begin_boundary
  acv%eend_boundary%previous%next=>acv%eend_boundary
 else
  acv%begin_boundary%next=>acv%eend_boundary
  acv%eend_boundary%previous=>acv%eend_boundary
 end if
 deallocate(acv%begin); deallocate(acv%eend)

 current=>acv%begin_boundary
 do while(associated(current))
  temp=>current%next
!   g_cell=current%g_cell
!   xxx=(associated(current))
  deallocate(current)
  current=>temp
 end do

end if

end subroutine delete_ac


subroutine delete_ac_node(temp,position)
! Delete a single node.

implicit none
type(auxi_crit_cell), pointer :: temp
! The pointer pointing to the will-be deleted critical cell.
character*6, intent(in) :: position
! Indicate whether the pointer will point to the critical cell
! previous or next to the deleted critical cell after the 
! delection.

type(auxi_crit_cell), pointer :: current
type(adss_info) :: adss

adss=temp%address
if(position.eq.'before') then
 if(associated(temp%previous)) then
  current=>temp%previous
 else
  current=>acvv(adss%cv_nmb)%begin
 end if
else
 if(associated(temp%next)) then
  current=>temp%next
 else
  current=>acvv(adss%cv_nmb)%eend
 end if
end if

if(associated(temp%previous)) then
 if(associated(temp%next)) then
  temp%previous%next=>temp%next
 else
  nullify(temp%previous%next)
 end if
else
 if(associated(temp%next)) then
  acvv(adss%cv_nmb)%begin%next=>temp%next
 else
  nullify(acvv(adss%cv_nmb)%begin%next)
 end if
end if
if(associated(temp%next)) then
 if(associated(temp%previous)) then
  temp%next%previous=>temp%previous
 else
  nullify(temp%next%previous)
 end if
else
 if(associated(temp%previous)) then
  acvv(adss%cv_nmb)%eend%previous=>temp%previous
 else
  nullify(acvv(adss%cv_nmb)%eend%previous)
 end if
end if     

if(acvv(adss%cv_nmb)%cv_type.eq.'circular') then
 if(acvv(adss%cv_nmb)%begin%next%address.eq.temp%address) then
  acvv(adss%cv_nmb)%begin%next=>temp%next
 end if
 if(acvv(adss%cv_nmb)%eend%previous%address.eq.temp%address) then
  acvv(adss%cv_nmb)%eend%previous=>temp%previous
 end if
end if

deallocate(temp)
temp=>current

end subroutine delete_ac_node


subroutine insert_ac(temp,temp_new,position)

implicit none
type(auxi_crit_cell), pointer :: temp, temp_new
character*6, intent(in) :: position

type(adss_info) :: ads

ads=temp%address
nullify(temp_new%previous); nullify(temp_new%next)
acvv(ads%cv_nmb)%total=acvv(ads%cv_nmb)%total+1
temp_new%address%idx=acvv(ads%cv_nmb)%total
select case(position)
 case('before')
  if(associated(temp%previous)) then
! The insert point is in the middle of the discontinuity curve.   
   temp_new%previous=>temp%previous
   temp_new%next=>temp
   temp%previous%next=>temp_new
   temp%previous=>temp_new
  else
   if(associated(acvv(ads%cv_nmb)%begin%next)) then
! The insert point is in front of the first critical cell in 
! the curve. 
    temp_new%next=>temp
	temp%previous=>temp_new
	acvv(ads%cv_nmb)%begin%next=>temp_new
   else
! The curve is still empty and 'temp_new' is going to be the 
! first member of the curve. 'temp' is the curve's 'begin'.
    acvv(ads%cv_nmb)%begin%next=>temp_new
	acvv(ads%cv_nmb)%eend%previous=>temp_new
   end if
  end if
 case('after ')
! The insert point is in front of the first critical cell in 
! the curve. 
  if(associated(temp%next)) then
   temp_new%previous=>temp
   temp_new%next=>temp%next
   temp%next%previous=>temp_new
   temp%next=>temp_new
  else
   if(associated(acvv(ads%cv_nmb)%eend%previous)) then
! The insert point is behind the last critical cell in the curve. 
    temp_new%previous=>temp
	temp%next=>temp_new
	acvv(ads%cv_nmb)%eend%previous=>temp_new
   else
! The curve is still empty and 'temp_new' is going ot be the 
! first member of the curve. 'temp' is the curve's 'eend'.
    acvv(ads%cv_nmb)%begin%next=>temp_new
	acvv(ads%cv_nmb)%eend%previous=>temp_new
   end if
  end if
end select

if(acvv(ads%cv_nmb)%cv_type.eq.'circular') then
 if(acvv(ads%cv_nmb)%begin%next%address.eq.temp%address) then
  if(position.eq.'before') acvv(ads%cv_nmb)%begin%next=>temp_new
 end if
 if(acvv(ads%cv_nmb)%eend%previous%address.eq.temp%address) then
  if(position.eq.'after') acvv(ads%cv_nmb)%eend%previous=>temp_new
 end if
end if

end subroutine insert_ac


subroutine creat_acv(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

allocate(acvv(cv_nmb)%begin); nullify(acvv(cv_nmb)%begin%previous)
nullify(acvv(cv_nmb)%begin%next); call clean_up_auxi_crit(acvv(cv_nmb)%begin)
allocate(acvv(cv_nmb)%eend); nullify(acvv(cv_nmb)%eend%next)
nullify(acvv(cv_nmb)%eend%previous)
call clean_up_auxi_crit(acvv(cv_nmb)%eend)
acvv(cv_nmb)%begin%address%cv_nmb=cv_nmb
acvv(cv_nmb)%eend%address%cv_nmb=cv_nmb

allocate(acvv(cv_nmb)%begin_boundary) 
nullify(acvv(cv_nmb)%begin_boundary%previous)
nullify(acvv(cv_nmb)%begin_boundary%next) 
call clean_up_auxi_crit(acvv(cv_nmb)%begin_boundary)
allocate(acvv(cv_nmb)%eend_boundary) 
nullify(acvv(cv_nmb)%eend_boundary%next)
nullify(acvv(cv_nmb)%eend_boundary%previous)
call clean_up_auxi_crit(acvv(cv_nmb)%eend_boundary)
acvv(cv_nmb)%begin_boundary%address%cv_nmb=cv_nmb
acvv(cv_nmb)%eend_boundary%address%cv_nmb=cv_nmb

end subroutine creat_acv


subroutine assign_curve(curve_2,curve_1)
! This subroutine defines the assignment of auxiliary discontinuity 
! curves.

implicit none
type(auxi_discv), intent(in) :: curve_1
type(auxi_discv), intent(out) :: curve_2

curve_2%status=curve_1%status
curve_2%cv_type=curve_1%cv_type
curve_2%begin_end=curve_1%begin_end
curve_2%end_end=curve_1%end_end
curve_2%total=curve_1%total
curve_2%wave=curve_1%wave
curve_2%begin=>curve_1%begin
curve_2%eend=>curve_1%eend

end subroutine assign_curve


subroutine side_area_length_partnered(g_cell,g_cell_p,l_area,r_area,length)

implicit none

type(geo_info), intent(in) :: g_cell, g_cell_p
real(8), intent(out) :: l_area, r_area, length

real(8), dimension(2,2) :: pts1, pts2
real(8), dimension(2) :: pt
character*10 :: err_msg
real(8) :: xxx

select case(g_cell%g_type)

 case('xx ')
  
  select case(g_cell_p%x_idx-g_cell%x_idx) 
   
   case(1) 
    
    if(g_cell%dis_posi(1).le.0.5d0.and.g_cell%dis_posi(2).le.0.5d0) then 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     xxx=0.5d0*(1.0d0-g_cell%dis_posi(1)-g_cell%dis_posi(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(g_cell%point(2).eq.'left  ') then
      l_area=xxx; r_area=1.0d0-xxx
     else
      r_area=xxx; l_area=1.0d0-xxx
     end if
     length=0.0d0
    else
     if(g_cell%dis_posi(1).gt.0.5d0.and.g_cell%dis_posi(2).gt.0.5d0) then 
      if(g_cell%point(2).eq.'left  ') then
       l_area=0.0d0; r_area=1.0d0
      else
       l_area=1.0d0; r_area=0.0d0
      end if
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      length=length_of_segment(pts1(1,:),pts1(2,:))
     else
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      pts2(:,1)=0.5d0; pts2(1,2)=-0.5d0; pts2(2,2)=0.5d0
      call ints_of_lines(pts1,pts2,pt,err_msg)
      if(g_cell%dis_posi(1).lt.0.5d0) then
       select case(g_cell%edge(1))
        case(1); xxx=0.5d0*(pt(2)+0.5d0)*(0.5d0-g_cell%dis_posi(1))
        case(3); xxx=0.5d0*(0.5d0-pt(2))*(0.5d0-g_cell%dis_posi(1))
        case default; call error_message
       end select
       length=length_of_segment(pts1(1,:),pt)
      else
       select case(g_cell%edge(2))
        case(1); xxx=0.5d0*(pt(2)+0.5d0)*(0.5d0-g_cell%dis_posi(2))
        case(3); xxx=0.5d0*(0.5d0-pt(2))*(0.5d0-g_cell%dis_posi(2))
        case default; call error_message
       end select
       length=length_of_segment(pts1(2,:),pt)
      end if
      if(g_cell%point(2).eq.'left  ') then
       l_area=xxx; r_area=1.0d0-xxx
      else
       l_area=1.0d0-xxx; r_area=xxx
      end if	  	  
     end if
    end if
    
   case(-1) 
    
    if(g_cell%dis_posi(1).ge.-0.5d0.and.g_cell%dis_posi(2).ge.-0.5d0) then 
     xxx=0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2)+1.0d0)
     if(g_cell%point(1).eq.'left  ') then
      l_area=xxx; r_area=1.0d0-xxx
     else
      r_area=xxx; l_area=1.0d0-xxx
     end if
    else
     if(g_cell%dis_posi(1).lt.-0.5d0.and.g_cell%dis_posi(2).lt.-0.5d0) then 
      if(g_cell%point(1).eq.'left  ') then
       l_area=0.0d0; r_area=1.0d0
      else
       l_area=1.0d0; r_area=0.0d0
      end if
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      length=length_of_segment(pts1(1,:),pts1(2,:))
     else
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      pts2(:,1)=-0.5d0; pts2(1,2)=-0.5d0; pts2(2,2)=0.5d0
      call ints_of_lines(pts1,pts2,pt,err_msg)
      if(g_cell%dis_posi(1).gt.-0.5d0) then
       select case(g_cell%edge(1))
        case(1); xxx=0.5d0*(pt(2)+0.5d0)*(g_cell%dis_posi(1)+0.5d0)
        case(3); xxx=0.5d0*(0.5d0-pt(2))*(g_cell%dis_posi(1)+0.5d0)
        case default; call error_message
       end select
       length=length_of_segment(pts1(1,:),pt)
      else
       select case(g_cell%edge(2))
        case(1); xxx=0.5d0*(pt(2)+0.5d0)*(g_cell%dis_posi(2)+0.5d0)
        case(3); xxx=0.5d0*(0.5d0-pt(2))*(g_cell%dis_posi(2)+0.5d0)
        case default; call error_message
       end select
       length=length_of_segment(pts1(2,:),pt)
      end if
      if(g_cell%point(1).eq.'left  ') then
       l_area=xxx; r_area=1.0d0-xxx
      else
       l_area=1.0d0-xxx; r_area=xxx
      end if	  	  
     end if
    end if
    
   case default; call error_message
  
  end select
  
 case('yy ')
  
  select case(g_cell_p%y_idx-g_cell%y_idx) 
   
   case(1) 
    
    if(g_cell%dis_posi(1).le.0.5d0.and.g_cell%dis_posi(2).le.0.5d0) then 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     xxx=0.5d0*(1.0d0-g_cell%dis_posi(1)-g_cell%dis_posi(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(g_cell%point(4).eq.'left  ') then
      l_area=xxx; r_area=1.0d0-xxx
     else
      r_area=xxx; l_area=1.0d0-xxx
     end if
    else
     if(g_cell%dis_posi(1).gt.0.5d0.and.g_cell%dis_posi(2).gt.0.5d0) then 
      if(g_cell%point(3).eq.'left  ') then
       l_area=0.0d0; r_area=1.0d0
      else
       l_area=1.0d0; r_area=0.0d0
      end if
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      length=length_of_segment(pts1(1,:),pts1(2,:))
     else
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      pts2(:,2)=0.5d0; pts2(1,1)=-0.5d0; pts2(2,1)=0.5d0
      call ints_of_lines(pts1,pts2,pt,err_msg)
      if(g_cell%dis_posi(1).lt.0.5d0) then
       select case(g_cell%edge(1))
        case(4); xxx=0.5d0*(pt(1)+0.5d0)*(0.5d0-g_cell%dis_posi(1))
        case(2); xxx=0.5d0*(0.5d0-pt(1))*(0.5d0-g_cell%dis_posi(1))
        case default; call error_message
       end select
       length=length_of_segment(pts1(1,:),pt)
      else
       select case(g_cell%edge(2))
        case(4); xxx=0.5d0*(pt(1)+0.5d0)*(0.5d0-g_cell%dis_posi(2))
        case(2); xxx=0.5d0*(0.5d0-pt(1))*(0.5d0-g_cell%dis_posi(2))
        case default; call error_message
       end select
       length=length_of_segment(pts1(2,:),pt)
      end if
      if(g_cell%point(4).eq.'left  ') then
       l_area=xxx; r_area=1.0d0-xxx
      else
       l_area=1.0d0-xxx; r_area=xxx
      end if	  	  
     end if
    end if
    
   case(-1) 
    
    if(g_cell%dis_posi(1).ge.-0.5d0.and.g_cell%dis_posi(2).ge.-0.5d0) then 
     xxx=0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2)+1.0d0)
     if(g_cell%point(1).eq.'left  ') then
      l_area=xxx; r_area=1.0d0-xxx
     else
      r_area=xxx; l_area=1.0d0-xxx
     end if
    else
     if(g_cell%dis_posi(1).lt.-0.5d0.and.g_cell%dis_posi(2).lt.-0.5d0) then 
      if(g_cell%point(1).eq.'left  ') then
       l_area=0.0d0; r_area=1.0d0
      else
       l_area=1.0d0; r_area=0.0d0
      end if
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      length=length_of_segment(pts1(1,:),pts1(2,:)) 
     else
      call pick_point(g_cell,1,pts1(1,:)); call pick_point(g_cell,2,pts1(2,:))
      pts2(:,2)=-0.5d0; pts2(1,1)=-0.5d0; pts2(2,1)=0.5d0
      call ints_of_lines(pts1,pts2,pt,err_msg)
      if(g_cell%dis_posi(1).gt.-0.5d0) then
       select case(g_cell%edge(1))
        case(4); xxx=0.5d0*(pt(1)+0.5d0)*(g_cell%dis_posi(1)+0.5d0)
        case(2); xxx=0.5d0*(0.5d0-pt(1))*(g_cell%dis_posi(1)+0.5d0)
        case default; call error_message
       end select
       length=length_of_segment(pts1(1,:),pt)
      else
       select case(g_cell%edge(2))
        case(4); xxx=0.5d0*(pt(1)+0.5d0)*(g_cell%dis_posi(2)+0.5d0)
        case(2); xxx=0.5d0*(0.5d0-pt(1))*(g_cell%dis_posi(2)+0.5d0)
        case default; call error_message
       end select
       length=length_of_segment(pts1(2,:),pt)
      end if
      if(g_cell%point(1).eq.'left  ') then
       l_area=xxx; r_area=1.0d0-xxx
      else
       l_area=1.0d0-xxx; r_area=xxx
      end if	  	  
     end if
    end if
    
   case default; call error_message
  
  end select
  
 case default; call error_message

end select

end subroutine side_area_length_partnered


end module auxiliary_discontinuity_curves