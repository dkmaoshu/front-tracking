module pick_middle_positions

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none

type(auxi_crit_cell), pointer :: tempa
type(geo_info) :: g_cell

public  pick_positions, tempa 
private g_cell


contains


subroutine pick_positions(up_down,positions,point_number, &
                          fashion,up_down_stream)

implicit none
character*6, intent(in) :: up_down
! The picking direction.
real(8), dimension(2,2), intent(out) :: positions
! The picked the discontinuity positions.
integer, intent(out) :: point_number
! The number of the points being picked.
character*8, intent(out) :: fashion
! The fashion of the picked points.
character*12, intent(out) :: up_down_stream
! The direction of the picked points.

type(auxi_crit_cell), pointer :: tempa_current
type(geo_info) :: g_current
real(8), dimension(2) :: pt, pt1, pt2
integer :: picking_number, i, num
! The number of the discontinuity position going to be picked in
! 'xy'-type critical cells.
logical :: if_length_ok, if_first_crossed, if_second_crossed, if_corner

g_cell=tempa%g_cell
positions=error_data
if_first_crossed=.false.; if_second_crossed=.false.
up_down_stream='no_given_yet'

! Deal with critical cells without previuos and next next-door
! neighbors.
if(.not.associated(tempa%p_nxt_dr).and..not.associated(tempa &
   %n_nxt_dr)) then
 if(g_cell%g_type.eq.'xy ') then
   
  do i=1, 2
   call check_crossed(g_cell,i,if_first_crossed)
   if(if_first_crossed) then
    point_number=1; fashion='single  '; return
   end if
  end do
  
! If 'crossed' happened to either middle discontinuity position, 
! the 'xy'-type critical cell should be handled forcedly.
  
  call check_length(g_cell,if_length_ok)
  if(.not.if_length_ok) then
   point_number=1; fashion='single  '; return
! If the length of the segment is too short, the critical cell 
! should also be handled forcedly.
  else
   call pick_middle_point(g_cell,1,pt1)
   call pick_middle_point(g_cell,2,pt2)
   select case(up_down)
    case('down  ')
     positions(1,:)=pt1; positions(2,:)=pt2; up_down_stream='down_stream '
    case('up    ')
     positions(1,:)=pt2; positions(2,:)=pt1; up_down_stream='up_stream   '
   end select
   point_number=2; fashion='one_side'; return  	   	 
! Otherwise, pick the two middle positions as the points.
  end if   
  
 else
! First determine the number of the starting discontinuity position.
  pt1=g_cell%mdis_posi; pt2=g_cell%mdis_posi
  select case(up_down)
   case('down  '); num=2
   case('up    '); num=1
   case default; call error_message
  end select

!  select case(up_down)
!   case('down  '); num=1
!   case('up    '); num=2
!   case default; call error_message
!  end select

! Select the two positions.
  select case(g_cell%g_type)
   case('xx ')
    select case(g_cell%edge(num))
     case(1); pt1(2)=-0.5d0; pt2(2)=0.5d0
     case(3); pt1(2)=0.5d0; pt2(2)=-0.5d0
     case default; call error_message
    end select
   case('yy ')
    select case(g_cell%edge(num))
     case(2); pt1(1)=-0.5d0; pt2(1)=0.5d0
     case(4); pt1(1)=0.5d0; pt2(1)=-0.5d0
     case default; call error_message
    end select
   case default; call error_message
  end select
! Produce the two positions.
  positions(1,:)=pt1; positions(2,:)=pt2;    		 	 	        
  select case(up_down)
   case('down  ')
    up_down_stream='down_stream '
   case('up    ')
    up_down_stream='up_stream   '
  end select

!  select case(up_down)
!   case('down  ')
!    positions(1,:)=pt1; positions(2,:)=pt2; up_down_stream='down_stream '
!   case('up    ')
!    positions(1,:)=pt2; positions(2,:)=pt1; up_down_stream='up_stream   '
!  end select

  point_number=2; fashion='one_side'; return
! For single 'xx'-, or 'yy'-type critical cells, both the two positions are chosen as
! the middle discontinuity positions.
 end if

end if

! Deal with the case of corner, i.e. two neighbored 'xy'-type
! critical cells with the same grid point as their single points.
if(g_cell%g_type.eq.'xy ') then
 call clean_up_g_cell(g_current)
 select case(up_down)
  case('up    ')
   if(associated(tempa%p_nxt_dr)) then
    g_current=tempa%p_nxt_dr%g_cell
    if(g_current%g_type.eq.'xy ') then
     call check_corner(g_cell,g_current,if_corner)
     if(if_corner) then
      fashion='corner  '; point_number=1; return
     else
      call error_message
     end if
    end if
   end if
  case('down  ')
   if(associated(tempa%n_nxt_dr)) then
    g_current=tempa%n_nxt_dr%g_cell
    if(g_current%g_type.eq.'xy ') then
     call check_corner(g_cell,g_current,if_corner)
     if(if_corner) then
      fashion='corner  '; point_number=1; return
     else
      call error_message
     end if
    end if
   end if
 end select 
end if

! Pick the first position.
select case(g_cell%g_type)

 case('xx ','yy ')
  positions(1,:)=g_cell%mdis_posi

 case('xy ')
  select case(up_down)
   case('up    ')
    call check_crossed(g_cell,2,if_first_crossed)
    if(.not.if_first_crossed) then
     call pick_middle_point(g_cell,2,pt)
    else
     call pick_middle_point(g_cell,1,pt)
    end if
   case('down  ')
    call check_crossed(g_cell,1,if_first_crossed)
    if(.not.if_first_crossed) then
     call pick_middle_point(g_cell,1,pt)
    else
     call pick_middle_point(g_cell,2,pt)
    end if
  end select
  positions(1,:)=pt

end select
point_number=1

! The following determines the picking direction and so on.
select case(up_down)
 
 case('up    ')
  if(associated(tempa%p_nxt_dr)) then
   tempa_current=>tempa%p_nxt_dr
   picking_number=1; fashion='between '; up_down_stream='up_stream   '
  else
   if(associated(tempa%n_nxt_dr)) then
    tempa_current=>tempa%n_nxt_dr
    picking_number=2; fashion='one_side'; up_down_stream='down_stream '
   else
    call error_message
   end if
  end if
      
 case('down  ')
  if(associated(tempa%n_nxt_dr)) then
   tempa_current=>tempa%n_nxt_dr
   picking_number=2; fashion='between '; up_down_stream='down_stream '
  else
   if(associated(tempa%p_nxt_dr)) then
    tempa_current=>tempa%p_nxt_dr
    picking_number=1; fashion='one_side'; up_down_stream='up_stream   '
   else
    call error_message
   end if
  end if
  
 case default; call error_message
       
end select

! Pick the second position.
g_current=tempa_current%g_cell
point_number=point_number+1

select case(g_current%g_type)
 case('xx ','yy ')
  pt=g_current%mdis_posi
 case('xy ')
  call check_crossed(g_current,picking_number,if_second_crossed)
  if(.not.if_second_crossed) then
   call pick_middle_point(g_current,picking_number,pt)
  else
   call pick_middle_point(g_current,3-picking_number,pt)
  end if
 case default; call error_message
end select

pt(1)=pt(1)+dfloat(g_current%x_idx-g_cell%x_idx)
pt(2)=pt(2)+dfloat(g_current%y_idx-g_cell%y_idx)
positions(2,:)=pt

! Modify the first position when 'fashion' is 'one_side'.
if(fashion.eq.'one_side'.and.g_cell%g_type.eq.'xy ') then
 select case(up_down)
  case('up    '); call pick_middle_point(g_cell,1,pt)
  case('down  '); call pick_middle_point(g_cell,2,pt)
 end select
 positions(1,:)=pt
end if

if(if_first_crossed) fashion='one_side'


contains


subroutine check_crossed(g_cell,num,if_crossed)

implicit none
type(geo_info), intent(in) :: g_cell
integer, intent(in) :: num
logical, intent(out) :: if_crossed

integer :: single

if(g_cell%g_type.ne.'xy ') call error_message
if_crossed=.false.

call find_single(g_cell,single)

select case(single)
 
 case(1)
  if(g_cell%dis_posi(num).lt.-0.5d0) then
   if_crossed=.true.; return
  end if
 
 case(2)
  if(g_cell%edge(num).eq.1.and.g_cell%dis_posi(num).gt.0.5d0.or. &
     g_cell%edge(num).eq.2.and.g_cell%dis_posi(num).lt.-0.5d0) then
   if_crossed=.true.; return
  end if

 case(3)
  if(g_cell%dis_posi(num).gt.0.5d0) then
   if_crossed=.true.; return
  end if
 
 case(4)
  if(g_cell%edge(num).eq.4.and.g_cell%dis_posi(num).gt.0.5d0.or. &
     g_cell%edge(num).eq.3.and.g_cell%dis_posi(num).lt.-0.5d0) then
   if_crossed=.true.; return
  end if

end select

end subroutine check_crossed


subroutine check_length(g_cell,if_ok)

implicit none
type(geo_info), intent(in) :: g_cell
logical, intent(out) :: if_ok

real(8), dimension(2) :: pt1, pt2
real(8) :: length

if_ok=.true.

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
length=length_of_segment(pt1,pt2)

if(length.lt.0.2d0) if_ok=.false.

end subroutine check_length


subroutine check_corner(g_cell1,g_cell2,if_corner)
! This subroutine checks whether the the two critical cells form a corner.

implicit none
type(geo_info), intent(in) :: g_cell1, g_cell2
logical, intent(out) :: if_corner

integer :: single1, single2, i

if_corner=.false.
call find_single(g_cell1,single1)
call find_single(g_cell2,single2)

if(iabs(g_cell1%x_idx-g_cell2%x_idx)+iabs(g_cell1%y_idx-g_cell2%y_idx).ne.1) then
 call error_message
end if

if(g_cell1%x_idx.ne.g_cell2%x_idx) then
 i=point_shift(single1,'xx ')
else
 i=point_shift(single1,'yy ')
end if
if(i.eq.single2) if_corner=.true.

end subroutine check_corner

   
end subroutine pick_positions

 
end module pick_middle_positions