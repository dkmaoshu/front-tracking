module compute_growth_rate

use discontinuity_curves

implicit none

real(8) :: tip_spike, tip_bubble

real(8):: hal_amplitude

real(8):: pre_amplitude

real(8):: growth_rate




private  cvv

public  tip_spike, tip_bubble

public pre_amplitude,hal_amplitude, get_growth_rate, growth_rate


contains


subroutine get_growth_rate(hal_amplitude,growth_rate)

implicit none

real(8):: hal_amplitude, pre_amplitude, growth_rate

real(8):: tip_bubble, tip_spike

type(critical_cell),pointer:: temp

allocate(cvv(curves_number))


call give_cv

temp=>cvv(1)%eend%previous

tip_bubble=(temp%g_cell%y_idx+temp%g_cell%dis_posi(2))*1.0/80.0

do while(abs(temp%g_cell%x_idx-40)>1.0d-12)
 temp=>temp%previous
end do

tip_spike=(temp%g_cell%y_idx+temp%g_cell%dis_posi(2))*1.0/80.0


hal_amplitude=0.5*(tip_spike-tip_bubble)

if(dabs(pre_amplitude).lt.1.0d-12) then
 growth_rate=0.0
else 
 growth_rate=(hal_amplitude-pre_amplitude)/(r*h)
end if

print*, growth_rate

pre_amplitude=hal_amplitude

call get_cv

deallocate(cvv)

end subroutine get_growth_rate

end module compute_growth_rate

