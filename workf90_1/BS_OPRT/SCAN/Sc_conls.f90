module scan_conservation_on_list

use	discontinuity_curves

implicit none

public  scan_con_list


contains


subroutine scan_con_list(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(state) :: al, ar, difference
real(8) :: aal, aar
integer :: head_mark
logical :: head_switch

if(cvv(cv_nmb)%status.eq.'asleep') then
 print*, 'The curve is still not waked up yet.'; return
end if

nullify(temp)
temp=>cvv(cv_nmb)%begin%next
head_mark=temp%address%idx; head_switch=.true.

do while(associated(temp).and.head_switch)

 g_cell=temp%g_cell; p_cell=temp%p_cell

 aal=side_area(g_cell,'left  ')
 al=aal*p_cell%l_state
 aar=side_area(g_cell,'right ')
 ar=aar*p_cell%r_state
 difference=p_cell%or_state-al-ar

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if

end do

end subroutine scan_con_list


end module scan_conservation_on_list
