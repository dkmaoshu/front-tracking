module auxiliary_list_fix
! The produced auxiliary lists may contain isolated coupled `xx'- or 
! `yy'-type auxiliary critical cells, i.e. the `xx'- or `yy'-type 
! auxiliary critical cells neighbored by `xy'-type auxiliary critical 
! cells at both ends. This module is for fixing this problem by
! deleting these type of coupled auxiliary critical cells.

use	auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use auxi_cc_prod
! 'auxi_ccp.f90'

implicit none

public list_fixing


contains


subroutine list_fixing

implicit none

integer :: cvnb

do cvnb=1, cvn
 if(acvv(cvnb)%status.eq.'awake ') then
  call single_isolated_fixing(cvnb)

!   call check_list_ac(1,'down  ')

  call coupled_isolated_fixing(cvnb)
  call nullify_orig(cvnb)
 end if
end do


contains


subroutine single_isolated_fixing(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: tempa
type(geo_info) :: g_cell, g_cell_xy_1, g_cell_xy_2

integer :: head_mark, ahead_number, behind_number
logical :: head_switch 
character*8 :: direction

tempa=>acvv(cv_nmb)%begin%next
head_mark=tempa%address%idx; head_switch=.true.

do while(associated(tempa).and.head_switch)

 call clean_up_g_cell(g_cell)
 call clean_up_g_cell(g_cell_xy_1); call clean_up_g_cell(g_cell_xy_2)

! Find the single isolated 'xx'- or 'yy'-critical cell.
 if(tempa%partner.ne.'single ') goto 10
 if(.not.associated(tempa%previous)) goto 10
 if(.not.associated(tempa%next)) goto 10
 g_cell=tempa%g_cell
 if(g_cell%g_type.ne.'xx '.and.g_cell%g_type.ne.'yy ') goto 10
 g_cell_xy_1=tempa%previous%g_cell
 g_cell_xy_2=tempa%next%g_cell
 if(g_cell_xy_1%g_type.ne.'xy '.or.g_cell_xy_2%g_type.ne.'xy ') &
   goto 10
  
! Find the numbers of coupled crtitical cells ahead and behind.
 call find_number(tempa,'ahead   ',ahead_number)
 call find_number(tempa,'behind  ',behind_number)
 direction='no_move '
 if(ahead_number.eq.2) then
  direction='ahead   '
 else
  if(behind_number.eq.2) then
   direction='behind  '
  else
   if(ahead_number.gt.2) then
    direction='ahead   '
   else
    if(behind_number.gt.2) direction='behind  '
   end if
  end if
 end if

! Fix the single isolated critical cell.
 if(direction.ne.'no_move ') call fixing(tempa,direction)
 
10  tempa=>tempa%next
 if(associated(tempa)) then
  head_switch=(tempa%address%idx.ne.head_mark)
 else
  exit
 end if

end do

end subroutine single_isolated_fixing


subroutine find_number(tempa,direction,number)

implicit none
type(auxi_crit_cell), pointer :: tempa
character*8, intent(in) :: direction
integer, intent(out) :: number

type(auxi_crit_cell), pointer :: current
type(geo_info) :: g_cell, g_cell_n

number=0
if(direction.eq.'ahead   ') then
  g_cell=tempa%next%g_cell; current=>tempa%next%next
else
 g_cell=tempa%previous%g_cell; current=>tempa%previous%previous
end if
do while(number.le.3)
 if(.not.associated(current).or.current%partner.eq.'single  ') exit
 g_cell_n=current%g_cell
 if(direction.eq.'ahead   ') then
  if(g_cell%edge(1).eq.g_cell_n%edge(2)) exit
 else
  if(g_cell%edge(2).eq.g_cell_n%edge(1)) exit
 end if   
 number=number+1
 if(direction.eq.'ahead   ') then
  g_cell=current%next%g_cell; current=>current%next%next
 else
  g_cell=current%previous%g_cell; current=>current%previous%previous
 end if
end do

end subroutine find_number


subroutine fixing(tempa,direction)

implicit none
type(auxi_crit_cell), pointer :: tempa
character*8, intent(in) :: direction

type(auxi_crit_cell), pointer :: current
type(geo_info) :: g_cell, g_cell_n, g_cell_nn
type(phy_info) :: p_cell, p_cell_n, p_cell_nn
type(geo_info), dimension(3) :: g_cell_fl
type(phy_info), dimension(3) :: p_cell_fl
integer :: i

call clean_up_g_cell(g_cell); call clean_up_p_cell(p_cell)
call clean_up_g_cell(g_cell_n); call clean_up_p_cell(p_cell_n)
call clean_up_g_cell(g_cell_nn); call clean_up_p_cell(p_cell_nn)
do i=1,3
 call clean_up_g_cell(g_cell_fl(i)); call clean_up_p_cell(p_cell_fl(i))
end do

select case(direction)
 case('ahead   ')
  current=>tempa%next
  g_cell=current%g_cell; p_cell=current%p_cell
  g_cell_n=current%next%temp%g_cell
  p_cell_n=current%next%temp%p_cell
  g_cell_nn=current%next%next%temp%g_cell
  p_cell_nn=current%next%next%temp%p_cell

  call with_geo(g_cell,g_cell_n,g_cell_fl(1),'next    ')
  call with_phy(g_cell,g_cell_n,p_cell,p_cell_n,p_cell_fl(1))
  call with_geo(g_cell_n,g_cell,g_cell_fl(2),'previous')
  call with_phy(g_cell_n,g_cell,p_cell_n,p_cell,p_cell_fl(2))
  g_cell_fl(3)=g_cell_nn; p_cell_fl(3)=p_cell_nn
 
  if(g_cell_fl(1)%g_type.eq.'non') return
  if(g_cell_fl(2)%g_type.eq.'non') call error_message

  current%g_cell=g_cell_fl(1); current%p_cell=p_cell_fl(1)
  current%partner='next    '

  current%next%g_cell=g_cell_fl(2)
  current%next%p_cell=p_cell_fl(2)
  current%next%partner='previous'

  current%next%next%g_cell=g_cell_fl(3)
  current%next%next%p_cell=p_cell_fl(3)
  current%next%next%partner='single  '

  current%p_nxt_dr=>tempa
  current%n_nxt_dr=>current%next%next
  
  current%next%p_nxt_dr=>tempa
  current%next%n_nxt_dr=>current%next%next

  current%next%next%p_nxt_dr=>current%next

  tempa=>current%next%next
 
 case('behind  ')
  current=>tempa%previous
  g_cell=current%g_cell; p_cell=current%p_cell
  g_cell_n=current%previous%temp%g_cell
  p_cell_n=current%previous%temp%p_cell
  g_cell_nn=current%previous%previous%temp%g_cell
  p_cell_nn=current%previous%previous%temp%p_cell

  call with_geo(g_cell_n,g_cell,g_cell_fl(2),'next    ')
  call with_phy(g_cell_n,g_cell,p_cell_n,p_cell,p_cell_fl(2))
  call with_geo(g_cell,g_cell_n,g_cell_fl(1),'previous')
  call with_phy(g_cell,g_cell_n,p_cell,p_cell_n,p_cell_fl(1))
  g_cell_fl(3)=g_cell_nn; p_cell_fl(3)=p_cell_nn

  if(g_cell_fl(1)%g_type.eq.'non') return
  if(g_cell_fl(2)%g_type.eq.'non') call error_message

  current%g_cell=g_cell_fl(1); current%p_cell=p_cell_fl(1)
  current%partner='previous'

  current%previous%g_cell=g_cell_fl(2)
  current%previous%p_cell=p_cell_fl(2)
  current%previous%partner='next'

  current%previous%previous%g_cell=g_cell_fl(3)
  current%previous%previous%p_cell=p_cell_fl(3)
  current%previous%previous%partner='single  '

  current%n_nxt_dr=>tempa
  current%p_nxt_dr=>current%previous%previous
  
  current%previous%n_nxt_dr=>tempa
  current%previous%p_nxt_dr=>current%previous%previous

  current%previous%previous%n_nxt_dr=>current%previous

 case default; call error_message

end select

end subroutine fixing


subroutine coupled_isolated_fixing(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: tempa
type(geo_info) :: g_cell, g_cell_n, g_cell_xy_1, g_cell_xy_2
type(phy_info) :: p_cell, p_cell_n, p_cell_xy_1, p_cell_xy_2
type(geo_info) :: g_cell_or, g_cell_orn
type(phy_info) :: p_cell_or, p_cell_orn
type(geo_info), dimension(4) :: g_cell_fl
type(phy_info), dimension(4) :: p_cell_fl

integer :: head_mark, i
logical :: head_switch !, if_isolated

tempa=>acvv(cv_nmb)%begin%next
head_mark=tempa%address%idx; head_switch=.true.

do while(associated(tempa).and.head_switch)

 call clean_up_g_cell(g_cell); call clean_up_g_cell(g_cell_n)
 call clean_up_g_cell(g_cell_xy_1); call clean_up_g_cell(g_cell_xy_2)
 call clean_up_p_cell(p_cell); call clean_up_p_cell(p_cell_n)
 call clean_up_p_cell(p_cell_xy_1); call clean_up_p_cell(p_cell_xy_2)
 call clean_up_g_cell(g_cell_or); call clean_up_p_cell(p_cell_or)
 call clean_up_g_cell(g_cell_orn); call clean_up_p_cell(p_cell_orn)
 do i=1,4
  call clean_up_g_cell(g_cell_fl(i)); call clean_up_p_cell(p_cell_fl(i))
 end do

! if_isolated=.true.

! First, find isolated auxiliary critical cell.
 if(tempa%partner.eq.'single  ') goto 10
 if(tempa%partner.eq.'previous') goto 10
 if(.not.associated(tempa%previous)) goto 10
 if(.not.associated(tempa%next%next)) goto 10
 g_cell=tempa%g_cell
 g_cell_n=tempa%next%g_cell
 g_cell_xy_1=tempa%previous%g_cell
 p_cell_xy_1=tempa%previous%p_cell
 g_cell_xy_2=tempa%next%next%g_cell
 p_cell_xy_2=tempa%next%next%p_cell
 if(g_cell_xy_1%g_type.ne.'xy '.or.g_cell_xy_2% &
    g_type.ne.'xy ') goto 10
 g_cell_or=tempa%temp%g_cell; p_cell_or=tempa%temp%p_cell
 g_cell_orn=tempa%next%temp%g_cell
 p_cell_orn=tempa%next%temp%p_cell
 if(g_cell_xy_1%edge(1).eq.g_cell_or%edge(2)) goto 10
 if(g_cell_xy_2%edge(2).eq.g_cell_orn%edge(1)) goto 10

! Fix the isolated auxiliary critical cell.
! if(if_isolated) then
  call with_geo(g_cell_xy_1,g_cell_or,g_cell_fl(1),'next    ')
  call with_phy(g_cell_xy_1,g_cell_or,p_cell_xy_1,p_cell_or, &
                p_cell_fl(1))
  call with_geo(g_cell_or,g_cell_xy_1,g_cell_fl(2),'previous')
  call with_phy(g_cell_or,g_cell_xy_1,p_cell_or,p_cell_xy_1, &
                p_cell_fl(2))
  call with_geo(g_cell_orn,g_cell_xy_2,g_cell_fl(3),'next    ')
  call with_phy(g_cell_orn,g_cell_xy_2,p_cell_orn,p_cell_xy_2, &
                p_cell_fl(3))
  call with_geo(g_cell_xy_2,g_cell_orn,g_cell_fl(4),'previous')
  call with_phy(g_cell_xy_2,g_cell_orn,p_cell_xy_2,p_cell_orn, &
                p_cell_fl(4))

  tempa%previous%g_cell=g_cell_fl(1)
  tempa%previous%p_cell=p_cell_fl(1)
  tempa%previous%partner='next    '  
  
  tempa%g_cell=g_cell_fl(2)
  tempa%p_cell=p_cell_fl(2)
  tempa%partner='previous'  
  
  tempa%next%g_cell=g_cell_fl(3)
  tempa%next%p_cell=p_cell_fl(3)
  tempa%next%partner='next    '  
  
  tempa%next%next%g_cell=g_cell_fl(4)
  tempa%next%next%p_cell=p_cell_fl(4)
  tempa%next%next%partner='previous'  

  tempa%previous%n_nxt_dr=>tempa%next
  tempa%p_nxt_dr=>tempa%previous%previous
  tempa%n_nxt_dr=>tempa%next
  tempa%next%p_nxt_dr=>tempa
  tempa%next%n_nxt_dr=>tempa%next%next%next
  tempa%next%next%p_nxt_dr=>tempa

  tempa=>tempa%next%next
! end if

10  tempa=>tempa%next
 if(associated(tempa)) then
  head_switch=(tempa%address%idx.ne.head_mark)
 else
  exit
 end if

end do


end subroutine coupled_isolated_fixing


subroutine nullify_orig(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: tempa
integer :: head_mark
logical :: head_switch

tempa=>acvv(cv_nmb)%begin%next
head_mark=tempa%address%idx; head_switch=.true.

do while(associated(tempa).and.head_switch)

 nullify(tempa%temp)

 tempa=>tempa%next
 if(associated(tempa)) then
  head_switch=(tempa%address%idx.ne.head_mark)
 else
  exit
 end if

end do

end subroutine nullify_orig


end subroutine list_fixing


end module auxiliary_list_fix

