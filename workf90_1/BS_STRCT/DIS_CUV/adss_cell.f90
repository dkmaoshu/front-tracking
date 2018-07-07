module adss_info_cell

implicit none

type adss_info
 integer :: cv_nmb
 integer :: idx
end type adss_info
! critical cell address information consists of the curve number 
! of the criticall, 'cv_nmd', and the index of the critical cell
! on the curve, 'idx'.

interface assignment(=)
 module procedure assignn
end interface

interface operator(.eq.)
 module procedure compare_true
end interface

interface operator(.ne.)
 module procedure compare_false
end interface

private compare_false, compare_true, assignn
public  clean_up_address


contains


subroutine assignn(address1,address2)

implicit none
type(adss_info), intent(out) :: address1
type(adss_info), intent(in) :: address2

address1%cv_nmb=address2%cv_nmb
address1%idx=address2%idx

end subroutine assignn


function compare_true(address1,address2) 

implicit none
type(adss_info), intent(in) :: address1, address2
logical :: compare_true

compare_true=.false.
if(address1%cv_nmb.eq.address2%cv_nmb) then
 if(address1%idx.eq.address2%idx) compare_true=.true.
end if

end function compare_true


function compare_false(address1,address2) 

implicit none
type(adss_info), intent(in) :: address1, address2
logical :: compare_false

compare_false=.false.
if(address1%cv_nmb.ne.address2%cv_nmb) then
 compare_false=.true.
else
 if(address1%idx.ne.address2%idx) compare_false=.true.
end if

end function compare_false


subroutine clean_up_address(address)

implicit none
type(adss_info), intent(inout) :: address

address%cv_nmb=0; address%idx=0

end subroutine clean_up_address

end module adss_info_cell