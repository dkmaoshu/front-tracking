module compress_linear_dis
! This module is for compressing linear discontinuites at the end 
! of a computation.

use solution
!'solution.f90'

implicit none

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(adss_info) :: address
integer :: k, steps, head_mark
logical :: head_switch

public  compress_linear_discontinuities
private temp, g_cell, p_cell, address, k, steps, head_mark, &
        head_switch


contains


subroutine compress_linear_discontinuities
! This subroutine compresses linear discontinuities.

implicit none

integer :: i

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ndd(ndn))

call give_sth(0)
call give_cv
call give_gcl
call give_nd

! call check_list_c(1,'down  ')

do i=1,cvn
 if(cvv(i)%status.eq.'awake ') then
  call produce_compressed_data(i)
  call update_with_compressed_data(i)

!   call check_list_c(1,'down  ')

 end if
end do

deallocate(uu); deallocate(cvv); deallocate(ggd_cell); deallocate(ndd)
 

contains


subroutine produce_compressed_data(cv_nb)
! This subroutine produces the compressed data.

implicit none

integer, intent(in) :: cv_nb

type(state) :: l_store, r_store

if(associated(cvv(cv_nb)%begin%next)) then
 temp=>cvv(cv_nb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch)
 g_cell=temp%g_cell; p_cell=temp%p_cell; address=temp%address
 if(p_cell%wv_nb.eq.2) then
  l_store=error_data; r_store=error_data
  if(temp%l_stk%cv_nmb.eq.0) then
   call pro_com_data('left  ',l_store)
  else
   l_store=p_cell%l_state
  end if
  if(temp%r_stk%cv_nmb.eq.0) then
   call pro_com_data('right ',r_store)
  else
   r_store=p_cell%r_state
  end if
! If the critical cell is stacked on a side, the compressing on the 
! side is not performed.  
  temp%l_store=l_store; temp%r_store=r_store
 end if

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine produce_compressed_data


subroutine pro_com_data(side,store)

implicit none

character*6, intent(in) :: side
type(state), intent(out) :: store

integer :: i0, j0, single, ii
type(state), dimension(-3:3) :: uxy
character*8 :: part
type(state) :: u_x, u_y

i0=g_cell%x_idx; j0=g_cell%y_idx
u_x=error_data; u_y=error_data

if(g_cell%g_type.eq.'xx'.or.g_cell%g_type.eq.'xy ') then
 call smth_dpr(i0,j0,temp,uxy,'xx ',side,'critical',.false.,3,'full    ')
 select case(g_cell%g_type)
  case('xx ')
   if(g_cell%point(1).eq.side) then
    part='negative' 
   else
    part='positive'
   end if    
  case('xy ')
   call find_single(g_cell,single)
   if(g_cell%point(single).eq.side) then
    select case(single)
     case(1,4); part='negative'
     case(2,3); part='positive'
    end select
   else
    select case(single)
     case(1,4); part='positive'
     case(2,3); part='negative'
    end select
   end if
 end select 	  
 select case(part)
  case('positive')
!   if(uxy(3).lt.0.9d0*error_data) then
   if(uxy(1).lt.0.9d0*error_data) then 
    u_x=uxy(0)
   else
    call inter_extrapolation(part,uxy,u_x)
   end if
  case('negative')
!   if(uxy(-3).lt.0.9d0*error_data) then
   if(uxy(-1).lt.0.9d0*error_data) then
    u_x=uxy(0)
   else
    call inter_extrapolation(part,uxy,u_x)
   end if
  case default; call error_message
 end select
end if 

if(g_cell%g_type.eq.'yy'.or.g_cell%g_type.ge.'xy ') then
 call smth_dpr(i0,j0,temp,uxy,'yy ',side,'critical',.false.,3,'full    ')
 select case(g_cell%g_type)
  case('yy ')
   if(g_cell%point(1).eq.side) then
    part='negative' 
   else
    part='positive'
   end if    
  case('xy ')
   call find_single(g_cell,single)
   if(g_cell%point(single).eq.side) then
    select case(single)
     case(1,2); part='negative'
     case(4,3); part='positive'
    end select
   else
    select case(single)
     case(1,2); part='positive'
     case(4,3); part='negative'
    end select
   end if
 end select 	  
 select case(part)
  case('positive')
!   if(uxy(3).lt.0.9d0*error_data) then
   if(uxy(1).lt.0.9d0*error_data) then 
    u_y=uxy(0)
   else
    call inter_extrapolation(part,uxy,u_y)
   end if
  case('negative')
!   if(uxy(-3).lt.0.9d0*error_data) then
   if(uxy(-1).lt.0.9d0*error_data) then 
    u_y=uxy(0)
   else
    call inter_extrapolation(part,uxy,u_y)
   end if
  case default; call error_message
 end select
end if 

if(u_x.gt.0.9d0*error_data.and.u_y.gt.0.9d0*error_data) then
 store=0.5d0*(u_x+u_y)
else
 if(u_x.gt.0.9d0*error_data) then
  store=u_x
 else
  if(u_y.gt.0.9d0*error_data) then
   store=u_y
  else
   call error_message
  end if
 end if
end if

end subroutine pro_com_data


subroutine inter_extrapolation(part,uxy,uu)

implicit none
character*8, intent(in) :: part
type(state), dimension(-3:3), intent(inout) :: uxy
type(state), intent(out) :: uu

integer :: stencil, i
logical :: if_ok
real(8) :: psu
type(state) :: uuu

select case(part)

 case('positive')
  if(uxy(1).gt.0.9d0*error_data) then
   if(uxy(2).gt.0.9d0*error_data) then
    if(uxy(3).gt.0.9d0*error_data) then
     uu=3.0d0*uxy(1)-3.0d0*uxy(2)+uxy(3)
     stencil=3
    else
     uu=2.0d0*uxy(1)-uxy(2)
     stencil=2
    end if
   else
    if(uxy(3).lt.0.9d0*error_data) then
     uu=uxy(1)
     stencil=1
    else
     call error_message
    end if
   end if
  else
   call error_message
  end if
  uuu=uxy(1) 
  
 case('negative')   
  if(uxy(-1).gt.0.9d0*error_data) then
   if(uxy(-2).gt.0.9d0*error_data) then
    if(uxy(-3).gt.0.9d0*error_data) then
     uu=3.0d0*uxy(-1)-3.0d0*uxy(-2)+uxy(-3)
     stencil=3
    else
     uu=2.0d0*uxy(-1)-uxy(-2)
     stencil=2
    end if
   else
    if(uxy(-1).gt.0.9d0*error_data) then
     uu=uxy(-1)
     stencil=1
    else
     call error_message
    end if
   end if
  else
   call error_message
  end if
  uuu=uxy(-1)
  
 case default; call error_message
end select

do i=1,4
 if(uxy(0)%value(i).gt.uuu%value(i)) then
  if(uu%value(i).gt.uxy(0)%value(i)) uu%value(i)=uxy(0)%value(i)
  if(uu%value(i).lt.uuu%value(i))  uu%value(i)=uuu%value(i)
 else
  if(uu%value(i).lt.uxy(0)%value(i)) uu%value(i)=uxy(0)%value(i)
  if(uu%value(i).gt.uuu%value(i))  uu%value(i)=uuu%value(i)
 end if
end do   

 if(dabs(uu%value(1)-uxy(0)%value(1)).lt.1.0d-14) then
  select case(part)
   case('positive'); uu%value(1)=0.9d0*uu%value(1)+0.1d0*uxy(1)%value(1)
   case('negative'); uu%value(1)=0.9d0*uu%value(1)+0.1d0*uxy(-1)%value(1)
   case default; call error_message
  end select
 end if

end subroutine inter_extrapolation


subroutine update_with_compressed_data(cv_nb)
! This subroutine updates the left and right states with the 
! compressed ones.

implicit none

integer, intent(in) :: cv_nb

type(state) :: l_physics, r_physics, l_conservative, r_conservative
real(8) :: l_area, r_area

if(associated(cvv(cv_nb)%begin%next)) then
 temp=>cvv(cv_nb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch)
 g_cell=temp%g_cell; p_cell=temp%p_cell; address=temp%address
 if(p_cell%wv_nb.eq.2) then
  l_area=side_area(g_cell,'left  ')
  r_area=side_area(g_cell,'right ')
  call find_physical_state(p_cell%l_state,(/1.0d0,0.0d0/),l_physics)
  call find_physical_state(p_cell%r_state,(/1.0d0,0.0d0/),r_physics)
  l_physics%value(1)=l_area*l_physics%value(1)+(1.0d0-l_area)*temp%l_store%value(1) 
  r_physics%value(1)=r_area*r_physics%value(1)+(1.0d0-r_area)*temp%r_store%value(1)
  call find_conservative_state(l_physics,(/1.0d0,0.0d0/),l_conservative)
  call find_conservative_state(r_physics,(/1.0d0,0.0d0/),r_conservative)
  p_cell%l_state=l_conservative
  p_cell%r_state=r_conservative
 end if

 temp%p_cell=p_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine update_with_compressed_data


end subroutine compress_linear_discontinuities


end module compress_linear_dis

