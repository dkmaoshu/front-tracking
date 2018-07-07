module produce_local_smooth_4_nodes
! This module produces local smooth solution for each two neighbor
! discontinuity curves pluging in the same node cell.

use produce_local_smooth
! 'pro_lcl.f90'

implicit none


private temp, temq, temp_a, temq_a, temp_p, temq_p, plug, plug_1, &
        prlcsm, i0, j0, side, side_1, ulc
public  prlcsmnd


contains


subroutine prlcsmnd(nd_nmb)
! This subroutine performs the production for the node with number
! 'nd_nmb'.

implicit none
integer, intent(in) :: nd_nmb

type(curve_plug), pointer :: temp_pp
integer :: head_mark, ii, jj
logical :: head_switch

i0=ndd(nd_nmb)%n_cell%x_idx; j0=ndd(nd_nmb)%n_cell%y_idx
temp_p=>ndd(nd_nmb)%cv_rg%begin%next
head_switch=.true.; head_mark=temp_p%address%idx
do while(head_switch)
! Produce local smooth solution in front of the current curve.
 call prlcsm(nd_nmb)

 do ii=-1,1; do jj=-1,1
  temp_p%plug%u_front(ii,jj)=ulc(ii,jj)
 end do; end do

 temp_p=>temp_p%next
 head_switch=(temp_p%address%idx.ne.head_mark)
end do

temp_p=>ndd(nd_nmb)%cv_rg%begin%next
head_switch=.true.; head_mark=temp_p%address%idx
do while(head_switch)
! Produce local smooth solution behind the current curve.
 
 temp_pp=>temp_p%previous

 do ii=-1,1; do jj=-1,1
  ulc(ii,jj)=temp_pp%plug%u_front(ii,jj)
  temp_p%plug%u_behind(ii,jj)=ulc(ii,jj)
 end do; end do

 temp_p=>temp_p%next
 head_switch=(temp_p%address%idx.ne.head_mark)
end do

end subroutine prlcsmnd 


subroutine prlcsm(nd_nmb)
! This subroutine produces local smooth solution in front of the
! current curve.

implicit none
integer, intent(in) :: nd_nmb
! The number of the node cell under concern.

! Initialization
nullify(temp); nullify(temp_a)
nullify(temq); nullify(temq_a)
nullify(temq_p)

side=temp_p%plug%side_front; plug=temp_p%plug
! Locate the head critical cell for the current discontinuity 
! curve 
select case(plug%end_type)

 case('begin ')
  select case(ndd(nd_nmb)%cv_rg%curves_type)
   case('regula')
    if(associated(cvv(plug%cv_nmb)%begin%next)) then
     temp=>cvv(plug%cv_nmb)%begin%next 
    end if
   case('auxila')
    if(associated(acvv(plug%cv_nmb)%begin%next)) then
     temp_a=>acvv(plug%cv_nmb)%begin%next 
    end if
  end select

 case('end   ')
  select case(ndd(nd_nmb)%cv_rg%curves_type)
   case('regula')
    if(associated(cvv(plug%cv_nmb)%eend%previous)) then
     temp=>cvv(plug%cv_nmb)%eend%previous 
    end if
   case('auxila')
    if(associated(acvv(plug%cv_nmb)%eend%previous)) then
     temp_a=>acvv(plug%cv_nmb)%eend%previous 
    end if
  end select
end select

! Find the next discontinuity curve pluging in the same node cell.
temq_p=>temp_p%next; plug_1=temq_p%plug

side_1=temq_p%plug%side_behind
! Locate the head critical cell for the next discontinuity curve
! pluging in the same node cell.
select case(plug_1%end_type)

 case('begin ')
  select case(ndd(nd_nmb)%cv_rg%curves_type)
   case('regula')
    if(associated(cvv(plug_1%cv_nmb)%begin%next)) then
     temq=>cvv(plug_1%cv_nmb)%begin%next 
    end if
   case('auxila')
    if(associated(acvv(plug_1%cv_nmb)%begin%next)) then
     temq_a=>acvv(plug_1%cv_nmb)%begin%next 
    end if
  end select

 case('end   ')
  select case(ndd(nd_nmb)%cv_rg%curves_type)
   case('regula')
    if(associated(cvv(plug_1%cv_nmb)%eend%previous)) then
     temq=>cvv(plug_1%cv_nmb)%eend%previous 
    end if
   case('auxila')
    if(associated(acvv(plug_1%cv_nmb)%eend%previous)) then
     temq_a=>acvv(plug_1%cv_nmb)%eend%previous 
    end if
  end select
end select

! Produce local smooth solution for each discontinuity curve 
! pluging in the node cell. The local smooth solution is in
! front ofated(tem each discontinuity curve.
select case(ndd(nd_nmb)%cv_rg%curves_type)
 case('regula'); call produce_smooth_regu
 case('auxila'); call produce_smooth_auxi
 case default; call error_message
end select

end subroutine prlcsm


end module produce_local_smooth_4_nodes