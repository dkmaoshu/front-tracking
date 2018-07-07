module local_smooth_fix
! Between two 'yawn'-status plugs the local smooth data have not been produced yet. This 
! module is for fixing this problem.

use solution
! 'solution.f90'

use solu_comput
! 'solu_com.f90'

implicit none

type(curve_plug), pointer :: tempp, tempp_n
type(cv_plug_info) :: plug, plug_n

private tempp, tempp_n, plug, plug_n
public  smnd_fix


contains


subroutine smnd_fix(nd_nmb)

implicit none
integer, intent(in) :: nd_nmb

type(curve_plug), pointer :: tempp_n1, tempp_n2
type(state), dimension(-1:1,-1:1) :: uu, uu1, uu2
integer :: head_mark, i, j, i0, j0, i1, j1, i2, j2
integer :: nd_nmb_n1, nd_nmb_n2
logical :: head_switch
character*6 :: status, status_n

do i=-1, 1; do j=-1, 1
 uu(i,j)=error_data
 uu1(i,j)=error_data 
 uu2(i,j)=error_data
end do; end do 
nullify(tempp); nullify(tempp_n); nullify(tempp_n1); nullify(tempp_n2)

i0=ndd(nd_nmb)%n_cell%x_idx; j0=ndd(nd_nmb)%n_cell%y_idx
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_switch=.true.; head_mark=tempp%address%idx
do while(head_switch)
 tempp_n=>tempp%next
 plug=tempp%plug; plug_n=tempp_n%plug

! Pick the status of the two neighboring discontinuity curves.
 select case(ndd(nd_nmb)%cv_rg%curves_type)
  case('regula')
   status=cvv(plug%cv_nmb)%status; status_n=cvv(plug_n%cv_nmb)%status
  case('auxila')
   status=acvv(plug%cv_nmb)%status; status_n=acvv(plug_n%cv_nmb)%status
  case default; call error_message
 end select

! It works only for the case that the two neighboring discontinuity curves are both
! 'yawn'. 
 if(status.eq.'yawn  '.and.status_n.eq.'yawn  ') then
  if(plug%u_front(0,0).gt.0.9*error_data) call error_message
  if(plug_n%u_behind(0,0).gt.0.9*error_data) call error_message

! Pick the two other node cells.  
  select case(ndd(nd_nmb)%cv_rg%curves_type)
   case('regula')
    if(plug%end_type.eq.'begin ') then
     nd_nmb_n1=cvv(plug%cv_nmb)%end_end
    else
     nd_nmb_n1=cvv(plug%cv_nmb)%begin_end
    end if
    if(plug_n%end_type.eq.'begin ') then
     nd_nmb_n2=cvv(plug_n%cv_nmb)%end_end
    else
     nd_nmb_n2=cvv(plug_n%cv_nmb)%begin_end
    end if
   case('auxila')
    if(plug%end_type.eq.'begin ') then
     nd_nmb_n1=acvv(plug%cv_nmb)%end_end
    else
     nd_nmb_n1=acvv(plug%cv_nmb)%begin_end
    end if
    if(plug_n%end_type.eq.'begin ') then
     nd_nmb_n2=acvv(plug_n%cv_nmb)%end_end
    else
     nd_nmb_n2=acvv(plug_n%cv_nmb)%begin_end
    end if
  end select
  i1=ndd(nd_nmb_n1)%n_cell%x_idx; j1=ndd(nd_nmb_n1)%n_cell%y_idx
  i2=ndd(nd_nmb_n2)%n_cell%x_idx; j2=ndd(nd_nmb_n2)%n_cell%y_idx

  if(iabs(i1-i0)+iabs(j1-j0).gt.1) call error_message
  if(iabs(i2-i0)+iabs(j2-j0).gt.1) call error_message

! Find the coresponding curve_plug of the curves on the othe node cells.
  if(plug%end_type.eq.'begin ') then
   call visit(nd_nmb_n1,plug%cv_nmb,'end   ',tempp_n1)
  else 
   call visit(nd_nmb_n1,plug%cv_nmb,'begin ',tempp_n1)
  end if
  do i=-1,1; do j=-1, 1
   uu1(i,j)=tempp_n1%plug%u_behind(i,j)
  end do; end do
  if(plug_n%end_type.eq.'begin ') then
   call visit(nd_nmb_n2,plug_n%cv_nmb,'end   ',tempp_n2)
  else 
   call visit(nd_nmb_n2,plug_n%cv_nmb,'begin ',tempp_n2)
  end if
  do i=-1, 1; do j=-1, 1
   uu2(i,j)=tempp_n2%plug%u_front(i,j)
  end do; end do

! Pick the smooth data.
  do i=-1, 1; do j=-1, 1
   if(i1-i0+i.le.1.and.i1-i0+i.ge.-1) then
    if(j1-j0+j.le.1.and.j1-j0+j.ge.-1) then
	 uu(i1-i0+i,j1-j0+j)=uu1(i,j)
    end if
   end if
  end do; end do
  do i=-1, 1; do j=-1, 1
   uu1(i,j)=uu(i,j); uu(i,j)=error_data
  end do; end do
  do i=-1, 1; do j=-1, 1
   if(i2-i0+i.le.1.and.i2-i0+i.ge.-1) then
    if(j2-j0+j.le.1.and.j2-j0+j.ge.-1) then
	 uu(i2-i0+i,j2-j0+j)=uu2(i,j)
    end if
   end if
  end do; end do
  do i=-1, 1; do j=-1, 1
   uu2(i,j)=uu(i,j); uu(i,j)=error_data
  end do; end do
     		  
! Finally, form the smooth data.
  do i=-1, 1; do j=-1, 1
   if(uu1(i,j).gt.0.9d0*error_data) then
    if(uu2(i,j).gt.0.9d0*error_data) then
     uu(i,j)=0.5d0*(uu1(i,j)+uu2(i,j))
    else
     uu(i,j)=uu1(i,j)
    end if
   else
    if(uu2(i,j).gt.0.9d0*error_data) then
     uu(i,j)=uu2(i,j)
    end if
   end if
  end do; end do   		 	   		  
    
  tempp%plug%u_front=uu; tempp_n%plug%u_behind=uu

 end if

 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)

end do

end subroutine smnd_fix

end module local_smooth_fix