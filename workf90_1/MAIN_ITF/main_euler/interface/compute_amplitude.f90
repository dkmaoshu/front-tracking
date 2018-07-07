Module compute_amplitude

use discontinuity_curves

implicit none

real(8) :: tip_spike, tip_bubble

public  tip_spike, tip_bubble, get_tip_of_spike, get_tip_of_bubble




contains

subroutine get_tip_of_spike
implicit none
type(critical_cell),pointer:: temp
temp=cvv(1)%eend
do while(temp%g_cell%x_idx>40)
 temp=>temp%previous
end do
tip_spike=temp%g_cell%dis_posi(2)
end subroutine get_tip_of_spike

subroutine get_tip_of_bubble
implicit none
type(critical_cell),pointer:: temp
temp=cvv(1)%eend
tip_bubble=temp%g_cell%dis_posi(2)
end subroutine get_tip_of_bubble

end module compute_amplitude
