module save_growth_rate

use compute_growth_rate

use input_setup

implicit none

real(8),dimension(1:2,1:100) :: amp
real(8),dimension(1:2,1:100) :: grt

public amp, grt

contains

subroutine save_gr(time_step,amp,grt)
implicit none

integer, intent(in) :: time_step

real(8),dimension(1:2,1:100) :: amp
real(8),dimension(1:2,1:100) :: grt

integer:: i

i=time_step
amp(1,i)=hal_amplitude
amp(2,i)=current_time
grt(1,i)=growth_rate
grt(2,i)=current_time

end subroutine save_gr
end module save_growth_rate
