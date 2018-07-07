module output_middle_positions
! This module is for output of middle discontinuity positions.
! The module is for tuning of the code.

use auxiliary_discontinuity_curves

implicit none

public  out_middle_posi


contains


subroutine out_middle_posi(lc)

implicit none
integer, intent(in) :: lc

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell
real(8), dimension(0:5000) :: x, y
integer :: head_mark, ll, i, l
logical :: head_switch

if(acvv(lc)%status.eq.'awake') then
 ll=0; temp=>acvv(lc)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 do while(associated(temp).and.head_switch)
  if(temp%partner.eq.'next    ') then
   temp=>temp%next
   cycle
  end if

  g_cell=temp%g_cell
  select case(g_cell%g_type)
   case('xx ','yy ')
    ll=ll+1
    x(ll)=(g_cell%x_idx+g_cell%mdis_posi(1))*h
    y(ll)=(g_cell%y_idx+g_cell%mdis_posi(2))*h
   case('xy ')
    do i=1,2
     ll=ll+1
     select case(g_cell%edge(i))
      case(1)
       x(ll)=(g_cell%x_idx+g_cell%mdis_posi(i))*h
       y(ll)=(g_cell%y_idx-0.5d0)*h
      case(2)
       x(ll)=(g_cell%x_idx+0.5d0)*h
       y(ll)=(g_cell%y_idx+g_cell%mdis_posi(i))*h
      case(3)
       x(ll)=(g_cell%x_idx+g_cell%mdis_posi(i))*h
       y(ll)=(g_cell%y_idx+0.5d0)*h
      case(4)
       x(ll)=(g_cell%x_idx-0.5d0)*h 
       y(ll)=(g_cell%y_idx+g_cell%mdis_posi(i))*h
     end select
    end do
  end select 
  temp=>temp%next
  if(associated(temp))head_switch=(temp%address%idx.ne.head_mark)
 end do

 open(8,file='d:\workf90_1\show\showmd.dat')
 do l=1,ll
  write(8,'(2f12.6)') x(l), y(l)
 end do
 close(8)

else
 
 print*, 'The curve does not exist yet.'
 pause

end if

end subroutine out_middle_posi


end module output_middle_positions
