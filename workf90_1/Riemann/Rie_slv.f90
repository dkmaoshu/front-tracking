program Riemann_solver

use euler_functions

implicit none

real(8), dimension(4) :: wb, wt, wmb, wmt
real(8) :: al, be, left_speed, right_speed

wb(1)=5.0d0; wb(2)=0.0d0; wb(3)=0.0d0; wb(4)=2.5d0
wt(1)=5.0d0; wt(2)=0.0d0; wt(3)=-23.6643d0; wt(4)=128.5d0

al=0.0d0; be=1.0d0

call triem(wb,wt,wmb,wmt,al,be)

left_speed=(flux_in_y(1,wb)-flux_in_y(1,wmb))/(wb(1)-wmb(1))
left_speed=(flux_in_y(3,wb)-flux_in_y(3,wmb))/(wb(3)-wmb(3))
left_speed=(flux_in_y(4,wb)-flux_in_y(4,wmb))/(wb(4)-wmb(4))

right_speed=(flux_in_y(1,wt)-flux_in_y(1,wmt))/(wt(1)-wmt(1))
right_speed=(flux_in_y(3,wt)-flux_in_y(3,wmt))/(wt(3)-wmt(3))
right_speed=(flux_in_y(4,wt)-flux_in_y(4,wmt))/(wt(4)-wmt(4))

end program Riemann_solver