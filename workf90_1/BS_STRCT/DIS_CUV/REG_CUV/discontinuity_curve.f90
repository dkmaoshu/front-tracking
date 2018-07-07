module discontinuity_curves

use grid, cvn=>curves_number
! 'grid.f90'

use critical_cells
! 'crit_cell.f90'

implicit none

type dis_curve
 type(critical_cell), pointer :: begin, eend
 character*6 :: status	 
! The possible values for variable 'status' are 'asleep',   
! 'awake' and 'yawn '.
 character*8 :: cv_type
! The possible values for variable 'cv_type' are 'dblink  ' and 
! 'circular'.
 integer :: begin_end, end_end
! The numbers of the nodes the two ends end; if zero, the end
! ends in infinitive boundary.
 integer :: total
! The total number of the critical cells on the curve.
 integer :: wave
end type dis_curve

type(dis_curve), dimension(cvn) :: cv
type(dis_curve), dimension(:), allocatable :: cvv

interface visit
 module procedure visit_c
end interface

interface deletee
 module procedure delete_c, delete_c_node
end interface

interface insert
 module procedure insert_c
end interface

private cv, delete_c, visit_c, delete_c_node, insert_c, &
        assign_curve
public  cvv, give_cv, get_cv, check_list_c, visit, deletee,  &
        insert, creat_cv, waking_curves, change_head_4_circular, &
        hidden_state_check_c

interface assignment(=)
 module procedure assign_curve
end interface


contains


subroutine give_cv

implicit none
integer :: i

! real(8), dimension(2) :: pt

do i=1,cvn
 cvv(i)%status=cv(i)%status
 cvv(i)%cv_type=cv(i)%cv_type
 cvv(i)%begin_end=cv(i)%begin_end; cvv(i)%end_end=cv(i)%end_end
 cvv(i)%total=cv(i)%total; cvv(i)%wave=cv(i)%wave
 nullify(cvv(i)%begin); nullify(cvv(i)%eend)
 if(cvv(i)%status.eq.'awake '.or.cvv(i)%status.eq.'yawn ') then
  cvv(i)%begin=>cv(i)%begin; cvv(i)%eend=>cv(i)%eend
 end if
end do

! call find_corner(pt,2)

end subroutine give_cv


subroutine get_cv

implicit none
integer :: i

do i=1,cvn
 cv(i)%status=cvv(i)%status
 cv(i)%cv_type=cvv(i)%cv_type
 cv(i)%begin_end=cvv(i)%begin_end; cv(i)%end_end=cvv(i)%end_end
 cv(i)%total=cvv(i)%total; cv(i)%wave=cvv(i)%wave
 nullify(cv(i)%begin); nullify(cv(i)%eend)
 if(cvv(i)%status.eq.'awake '.or.cvv(i)%status.eq.'yawn ') then
  cv(i)%begin=>cvv(i)%begin; cv(i)%eend=>cvv(i)%eend
 end if
end do

end subroutine get_cv


subroutine check_list_c(cv_nmb,up_down)

implicit none

character*6, intent(in) :: up_down
integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp, temp_g
type(geo_info) :: gc_c, gc_p, gc_n, gc_l, gc_r
type(phy_info) :: pc_c, pc_p, pc_n, pc_l, pc_r
type(phy_info) :: phy_state
type(adss_info) :: address
integer :: head_mark
logical :: head_switch
real(8), dimension(2) :: normal, pt1, pt2
  
! real(8) :: xxx, yyy, aaa
  
! xxx=0.0d0; yyy=0.0d0 
  
if(cvv(cv_nmb)%status.eq.'asleep') then
 print*, 'The curve is still not waked up yet.'; return
end if

nullify(temp)
if(up_down.eq.'up') then
 if(associated(cvv(cv_nmb)%eend%previous)) then
  temp=>cvv(cv_nmb)%eend%previous
  head_mark=temp%address%idx; head_switch=.true.
 else
  print*, 'The curve list is empty.'; return
 end if    
else
 if(associated(cvv(cv_nmb)%begin%next)) then
  temp=>cvv(cv_nmb)%begin%next
  head_mark=temp%address%idx; head_switch=.true.
 else
  print*, 'The curve list is empty.'; return
 end if
end if  

do while(associated(temp).and.head_switch)
 
 call clean_up_g_cell(gc_c); call clean_up_g_cell(gc_p); call clean_up_g_cell(gc_n)
 call clean_up_g_cell(gc_l); call clean_up_g_cell(gc_r)
 call clean_up_p_cell(pc_c); call clean_up_p_cell(pc_p); call clean_up_p_cell(pc_n)
 call clean_up_p_cell(pc_l); call clean_up_p_cell(pc_r)

 gc_c=temp%g_cell; pc_c=temp%p_cell; address=temp%address

 if(associated(temp%previous)) then
  gc_p=temp%previous%g_cell; pc_p=temp%previous%p_cell
 endif
 if(associated(temp%next)) then
  gc_n=temp%next%g_cell; pc_n=temp%next%p_cell
 end if
 if(temp%l_stk%cv_nmb.gt.0) then
  call visit(temp%l_stk,temp_g)
  gc_l=temp_g%g_cell; pc_l=temp_g%p_cell
 end if 
 if(temp%r_stk%cv_nmb.gt.0) then
  call visit(temp%r_stk,temp_g)

  gc_r=temp_g%g_cell; pc_r=temp_g%p_cell
 end if

 call pick_point(gc_c,1,pt1); call pick_point(gc_c,2,pt2)
 normal=normal_of_line(pt1,pt2)
  
 normal=(/1.0d0,0.0d0/)

 call find_physical_state(pc_c%l_state,normal,phy_state%l_state)
 call find_physical_state(pc_c%r_state,normal,phy_state%r_state)
! call find_physical_state(pc_c%or_state,normal,phy_state%or_state)
    
!   xxx=internal_energy(pc_c%or_state)
!   if(xxx.lt.0.0d0) call error_message
    
!   if(gc_c%g_type.eq.'xx ') then
!    if(gc_c%edge(1).eq.1) then
!     if(dfloat(gc_c%x_idx)*h.gt.1.0d0) then 
!      if(gc_c%dis_posi(1).lt.gc_c%dis_posi(2)) print*, gc_c%x_idx, gc_c%y_idx
!     else
!      if(gc_c%dis_posi(1).gt.gc_c%dis_posi(2)) print*, gc_c%x_idx, gc_c%y_idx
!     end if
!    else
!     if(dfloat(gc_c%x_idx)*h.gt.1.0d0) then 
!      if(gc_c%dis_posi(1).gt.gc_c%dis_posi(2)) print*, gc_c%x_idx, gc_c%y_idx
!     else
!      if(gc_c%dis_posi(1).lt.gc_c%dis_posi(2)) print*, gc_c%x_idx, gc_c%y_idx
!     end if
!    end if
!   end if
    
!   xxx=dfloat(gc_c%x_idx)*h-1.0d0+0.5d0*(gc_c%dis_posi(1)+gc_c%dis_posi(2))*h
!   
!   xxx=1.0d0-dfloat(gc_c%x_idx)*h-0.5d0*(gc_c%dis_posi(1)+gc_c%dis_posi(2))*h
   
!   xxx=(dfloat(gc_c%y_idx)+0.5d0*(gc_c%dis_posi(1)+gc_c%dis_posi(2)))*h-49.50d0*h
   
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Compute spanwise width for Air-R22.
!    
!  call pick_point(gc_c,1,pt1)
!  call pick_point(gc_c,2,pt2)
!  aaa=(dfloat(gc_c%y_idx)+dmax1(pt1(2),pt2(2)))*h
!  if(aaa.gt.xxx) then
!   xxx=aaa; print*,'+ ', gc_c%x_idx, gc_c%y_idx
!  end if
!  aaa=(dfloat(gc_c%y_idx)+dmin1(pt1(2),pt2(2)))*h
!  if(aaa.lt.yyy) then
!   yyy=aaa; print*,'- ', gc_c%x_idx, gc_c%y_idx
!  end if
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  
! continue
 
end subroutine check_list_c


subroutine visit_c(adss,temp_g)
! The subroutine is for visit of critical cells by their addresses.
! Variable 'adss' is the address of the critical cell to be 
! visited, and variable 'temp_g' is the guest pointer of
! critical cell containing the information passed in or out.

implicit none

type(adss_info), intent(in)	:: adss
type(critical_cell), pointer :: temp_g

type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: l_stk, r_stk


type(critical_cell), pointer :: temp
integer :: head_mark
logical :: head_switch

if(cvv(adss%cv_nmb)%status.eq.'asleep ') then
 print*, 'The curve list does not exist.'; return
end if

if(associated(cvv(adss%cv_nmb)%begin%next)) then
 temp=>cvv(adss%cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if

nullify(temp_g)

do while(associated(temp).and.head_switch) 

 if(temp%address%idx.eq.adss%idx) then
  temp_g=>temp
  g_cell=temp%g_cell; p_cell=temp%p_cell
  l_stk=temp%l_stk; r_stk=temp%r_stk
  exit
 endif       
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  print*, 'No match in the visit.'
  call error_message
 endif

end do   

if(.not.associated(temp_g)) then
 print*, 'No match in the visit.'
 call error_message
end if

end subroutine visit_c


subroutine hidden_state_check_c

implicit none

integer :: cv_nmb
type(critical_cell), pointer :: temp, temp_n
type(geo_info) :: g_cell, g_cell_n
type(phy_info) :: p_cell, p_cell_n 
integer :: head_mark
logical :: head_switch
real(8) :: difference

do cv_nmb=1,cvn
 if(cvv(cv_nmb)%status.eq.'asleep') return
  
 nullify(temp); nullify(temp_n)
  
 if(associated(cvv(cv_nmb)%begin%next)) then
  temp=>cvv(cv_nmb)%begin%next
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

end subroutine hidden_state_check_c


subroutine delete_c(cvv)

implicit none
type(dis_curve), intent(inout) :: cvv

type(critical_cell), pointer :: temp, current

if(cvv%status.eq.'awake '.or.cvv%status.eq.'yawn  ') then
 cvv%status='asleep'
 cvv%cv_type=' '
 cvv%begin_end=-1000; cvv%end_end=-1000
 cvv%total=-1000

 if(associated(cvv%begin%next)) then
  cvv%begin%next%previous=>cvv%begin
  cvv%eend%previous%next=>cvv%eend
 else
  deallocate(cvv%begin)
  deallocate(cvv%eend)
  return
 end if

 current=>cvv%begin
 do while(associated(current)) 
  temp=>current%next
  deallocate(current)
  current=>temp
 end do

end if

end subroutine delete_c


subroutine delete_c_node(temp,position)

implicit none
type(critical_cell), pointer :: temp
! The pointer pointing to the will-be deleted critical cell.
character*6, intent(in) :: position
! Indicate whether the pointer will point to the critical cell
! previous or next to the deleted critical cell after the
! delection.

type(critical_cell), pointer :: current
type(adss_info) :: adss

adss=temp%address
if(position.eq.'before') then
 if(associated(temp%previous)) then
  current=>temp%previous
 else
  current=>cvv(adss%cv_nmb)%begin
 end if
else
 if(associated(temp%next)) then
  current=>temp%next
 else
  current=>cvv(adss%cv_nmb)%eend
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
  cvv(adss%cv_nmb)%begin%next=>temp%next
 else
  nullify(cvv(adss%cv_nmb)%begin%next)
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
  cvv(adss%cv_nmb)%eend%previous=>temp%previous
 else
  nullify(cvv(adss%cv_nmb)%eend%previous)
 end if
end if     

if(cvv(adss%cv_nmb)%cv_type.eq.'circular') then
 if(cvv(adss%cv_nmb)%begin%next%address.eq.temp%address) then
  cvv(adss%cv_nmb)%begin%next=>temp%next
 end if
 if(cvv(adss%cv_nmb)%eend%previous%address.eq.temp%address) then
  cvv(adss%cv_nmb)%eend%previous=>temp%previous
 end if
end if

deallocate(temp)
temp=>current

end subroutine delete_c_node


subroutine insert_c(temp,temp_new,position)

implicit none
type(critical_cell), pointer :: temp, temp_new
character*6, intent(in) :: position

type(geo_info) :: g_cell, g_new
type(adss_info) :: ads, ads_new

ads=temp%address
g_cell=temp%g_cell; g_new=temp_new%g_cell
nullify(temp_new%previous); nullify(temp_new%next)
cvv(ads%cv_nmb)%total=cvv(ads%cv_nmb)%total+1
temp_new%address%cv_nmb=temp%address%cv_nmb
temp_new%address%idx=cvv(ads%cv_nmb)%total
temp_new%l_stk=adss_info(0,0); temp_new%r_stk=adss_info(0,0)
ads_new=temp_new%address
select case(position)
 case('before')
  if(associated(temp%previous)) then
! The insert point is in the middle of the discontinuity curve.   
   temp_new%previous=>temp%previous
   temp_new%next=>temp
   temp%previous%next=>temp_new
   temp%previous=>temp_new
  else
   if(associated(cvv(ads%cv_nmb)%begin%next)) then
! The insert point is in front of the first critical cell in 
! the curve. 
    temp_new%next=>temp
    temp%previous=>temp_new
    cvv(ads%cv_nmb)%begin%next=>temp_new
   else
! The curve is still empty and 'temp_new' is going to be the 
! first member of the curve. 'temp' is the curve's 'begin'.
    cvv(ads%cv_nmb)%begin%next=>temp_new
    cvv(ads%cv_nmb)%eend%previous=>temp_new
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
   if(associated(cvv(ads%cv_nmb)%eend%previous)) then
! The insert point is behind the last critical cell in the curve. 
    temp_new%previous=>temp
	temp%next=>temp_new
	cvv(ads%cv_nmb)%eend%previous=>temp_new
   else
! The curve is still empty and 'temp_new' is going ot be the 
! first member of the curve. 'temp' is the curve's 'eend'.
    cvv(ads%cv_nmb)%begin%next=>temp_new
	cvv(ads%cv_nmb)%eend%previous=>temp_new
   end if
  end if
end select

if(cvv(ads%cv_nmb)%cv_type.eq.'circular') then
 if(cvv(ads%cv_nmb)%begin%next%address.eq.temp%address) then
  if(position.eq.'before') cvv(ads%cv_nmb)%begin%next=>temp_new
 end if
 if(cvv(ads%cv_nmb)%eend%previous%address.eq.temp%address) then
  if(position.eq.'after ') cvv(ads%cv_nmb)%eend%previous=>temp_new
 end if
end if

end subroutine insert_c


subroutine creat_cv(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

allocate(cvv(cv_nmb)%begin); nullify(cvv(cv_nmb)%begin%previous)
nullify(cvv(cv_nmb)%begin%next); call clean_up_crit(cvv(cv_nmb)%begin)
allocate(cvv(cv_nmb)%eend); nullify(cvv(cv_nmb)%eend%next)
nullify(cvv(cv_nmb)%eend%previous); call clean_up_crit(cvv(cv_nmb)%eend)

cvv(cv_nmb)%begin%address%cv_nmb=cv_nmb
cvv(cv_nmb)%eend%address%cv_nmb=cv_nmb

end subroutine creat_cv


subroutine waking_curves
! This subrouitine wakes all the discontinuity curves in 'yawn  ' 
! mode.

implicit none
integer :: i

do i=1, cvn
 if(cvv(i)%status.eq.'yawn  ') then
  if(associated(cvv(i)%begin%next)) cvv(i)%status='awake '
 end if
end do

end subroutine waking_curves


subroutine assign_curve(curve_2,curve_1)
! This subroutine defines the assignment of discontinuity curves.

implicit none
type(dis_curve), intent(in) :: curve_1
type(dis_curve), intent(out) :: curve_2

curve_2%status=curve_1%status
curve_2%cv_type=curve_1%cv_type
curve_2%begin_end=curve_1%begin_end
curve_2%end_end=curve_1%end_end
curve_2%total=curve_1%total
curve_2%wave=curve_1%wave
curve_2%begin=>curve_1%begin
curve_2%eend=>curve_1%eend

end subroutine assign_curve


subroutine put_empty_yawn(cv_nmb)
! This subroutine puts empty discontinuity curve list with 'awake'
! status back to 'yawn' status.

implicit none
integer, intent(in) :: cv_nmb

if(cvv(cv_nmb)%status.eq.'awake ') then
 if(.not.associated(cvv(cv_nmb)%begin%next)) then
  if(.not.associated(cvv(cv_nmb)%eend%previous)) then
   cvv(cv_nmb)%status='yawn  '
  end if
 end if
end if
 
end subroutine put_empty_yawn


subroutine change_head_4_circular
! It may not be good that the head critical cell of a circular
! discontinuity curve is of 'xy'-type. In this case, we should 
! shift the head to an 'xx'- or 'yy'-type critical cell.

implicit none
type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
integer :: head_mark, l
logical :: head_switch
character*3 :: found

do l=1, cvn
 if(cvv(l)%status.eq.'awake '.and.cvv(l)% &
    cv_type.eq.'circular') then
  temp=>cvv(l)%begin%next; g_cell=temp%g_cell
  if(g_cell%g_type.eq.'xy ') then
   head_mark=temp%address%idx; head_switch=.true.
   found='no '
   do while(head_switch)
    temp=>temp%next; g_cell=temp%g_cell
    if(g_cell%g_type.ne.'xy ') then
     found='yes'; exit
    end if
    if(temp%address%idx.eq.head_mark) head_switch=.false.
   end do
   if(found.eq.'no ') call error_message
   cvv(l)%begin%next=>temp
   cvv(l)%eend%previous=>temp%next
  end if
 end if
end do

end subroutine change_head_4_circular


end module discontinuity_curves									   