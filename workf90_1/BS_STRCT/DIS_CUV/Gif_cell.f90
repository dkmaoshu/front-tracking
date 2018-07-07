module geo_info_cell
! This module describes the data structure of the geometrical 
! information contained in critical cells.

use grid
! 'grid.f90'

use tools
! 'tools.f90'

type geo_info

 character*3 :: g_type	   
 integer :: x_idx, y_idx   
 real(8) :: dis_posi(2)
 real(8) :: mdis_posi(2)
 real(8), dimension(2,2) :: normal
 integer :: edge(2)
 character*3 :: wind(2)
 character*6 :: point(4)
end type geo_info
! Geometrical information consists of the type of a critical cell
! 'c_type', the x- and y-indices of the critical cell, the two 
! discontinuity positions 'dis_posi', x- and y-coordinates of
! the middle discontinuity positions 'mdis_posi', the normals of 
! the curve at the two discontinuity positionstwo cell-edges 
! that crossed by the two discontinuity positions 'edge' and the
! wind direction on these two edges, and the sides the four corner 
! points belong 'point'. 

! The orientation of the discontinuity curve in the critical cell
! is defined as pointing from the first discontinuity position to
! the second one.

public find_points_from_edges, find_edges_from_points, &
find_type_from_edges, find_type_from_points, clean_up_g_cell, &
point_shift, check_smooth, pick_point, find_single, &
pick_corner, side_shift, side_area, compute_curve_length, &
pick_middle_point, neighbor_edge, which_side, normal_4_isolated_xy, &
side_of_point, length_of_curve_segment


contains


subroutine find_points_from_edges(g_cell)

implicit none
type(geo_info), intent(inout) :: g_cell

integer :: egn(3), ptn(4), i, j

do i=1,3;  egn(i)=cycle_4(g_cell%edge(1)+i); end do
do i=1,4;  ptn(i)=cycle_4(g_cell%edge(1)+i); end do

do i=1,3
 if(g_cell%edge(2).eq.egn(i)) then
  do j=1, i; g_cell%point(ptn(j))='right'; end do
  do j=i+1, 4; g_cell%point(ptn(j))='left'; end do
  exit
 end if
end do

end subroutine find_points_from_edges


subroutine find_edges_from_points(g_cell)

implicit none
type(geo_info), intent(inout) :: g_cell

integer :: i, j

do i=1,4
 j=cycle_4(i+1)
 if(g_cell%point(i).ne.g_cell%point(j)) then
  if(g_cell%point(i).eq.'left') then
   g_cell%edge(1)=i
  else
   g_cell%edge(2)=i
  end if
 end if      
end do

end subroutine find_edges_from_points


subroutine find_type_from_edges(g_cell)

implicit none
type(geo_info),intent(inout) :: g_cell

if(g_cell%edge(1).eq.1.or.g_cell%edge(1).eq.3) then
 if(g_cell%edge(2).eq.1.or.g_cell%edge(2).eq.3) then
  g_cell%g_type='xx'
 else
  g_cell%g_type='xy'
 end if
end if
if(g_cell%edge(1).eq.2.or.g_cell%edge(1).eq.4) then
 if(g_cell%edge(2).eq.2.or.g_cell%edge(2).eq.4) then
  g_cell%g_type='yy'
 else
  g_cell%g_type='xy'
 end if
end if

end subroutine find_type_from_edges


subroutine find_type_from_points(g_cell)

implicit none
type(geo_info),intent(inout) :: g_cell

call find_edges_from_points(g_cell)
call find_type_from_edges(g_cell)

end subroutine find_type_from_points


function side_shift(side)
implicit none
character*6, intent(in) :: side
character*6 :: side_shift
if(side.eq.'left') then
 side_shift='right'
else
 if(side.eq.'right') then
  side_shift='left'
 else
  print*, 'The input "side" is incorrectly definded'
  stop
 end if
end if
 
end function side_shift

function point_shift(point,direction)
implicit none        
integer, intent(in) :: point
character*3, intent(in) :: direction
integer :: point_shift

select case(point)
 case(1)
  select case(direction)
   case('xx'); point_shift=2
   case('yy'); point_shift=4
   case('xy'); point_shift=3
   case default; call error_message
  end select
 case(2)
  select case(direction)
   case('xx'); point_shift=1
   case('yy'); point_shift=3
   case('xy'); point_shift=4
   case default; call error_message
  end select
 case(3)
  select case(direction)
   case('xx'); point_shift=4
   case('yy'); point_shift=2
   case('xy'); point_shift=1
   case default; call error_message
  end select
 case(4)
  select case(direction)
   case('xx'); point_shift=3
   case('yy'); point_shift=1
   case('xy'); point_shift=2
   case default; call error_message
  end select
 case default; call error_message
end select

end function point_shift   	  	 


subroutine clean_up_g_cell(g_cell)
! Set default value for a geometry cell of the auxiliary 
! critical cell to be produced.

implicit none

type(geo_info), intent(out) :: g_cell

g_cell%g_type='non'
g_cell%x_idx=-1000; g_cell%y_idx=-1000
g_cell%dis_posi=error_data
g_cell%mdis_posi=error_data
g_cell%normal=error_data
g_cell%edge=-1000
g_cell%point=' '

end subroutine clean_up_g_cell


function neighbor_edge(edge)

implicit none
integer, intent(in) :: edge
integer :: neighbor_edge

select case(edge)
 case(1); neighbor_edge=3
 case(2); neighbor_edge=4
 case(3); neighbor_edge=1
 case(4); neighbor_edge=2
end select

end function neighbor_edge


subroutine check_smooth(g_cell,gn_cell,smooth)
! This subroutine checks if the discontinuity curve is smooth 
! near the given critical cells. The check is done on the 
! observation of the angle between the two normals of the curve
! segments in the two critical cell. If the angle is greater than
! '0.1*pi', it is not smooth.

implicit none

type(geo_info), intent(in) :: g_cell, gn_cell
character*3, intent(out) :: smooth

real(8), dimension(2) :: normal_1, normal_2
real(8) :: angle
real(8), dimension(2) :: pt1, pt2

smooth='yes'

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
normal_1=normal_of_line(pt1,pt2)

call pick_point(gn_cell,1,pt1)
call pick_point(gn_cell,2,pt2)
normal_2=normal_of_line(pt1,pt2)

angle=normal_1(1)*normal_2(2)-normal_2(1)*normal_1(2)
angle=dabs(angle)
angle=dasin(angle)

if(angle.gt.3.1415/6.0d0) smooth='no'

end subroutine check_smooth


subroutine pick_point(g_cell,num,pt)
! This subroutine picks the point corresponding to the 'num'th
! discontinuity position in a critical cell.

implicit none

type(geo_info), intent(in) :: g_cell
integer, intent(in) :: num
real(8), dimension(2), intent(out) :: pt

select case(g_cell%edge(num))
 case(1); pt(1)=g_cell%dis_posi(num); pt(2)=-0.5d0
 case(2); pt(1)=0.5d0; pt(2)=g_cell%dis_posi(num)
 case(3); pt(1)=g_cell%dis_posi(num); pt(2)=0.5d0
 case(4); pt(1)=-0.5d0; pt(2)=g_cell%dis_posi(num)
end select

end subroutine pick_point


subroutine pick_middle_point(g_cell,num,pt)
! This subroutine picks the point corresponding to the 'num'th
! discontinuity position in a critical cell.

implicit none

type(geo_info), intent(in) :: g_cell
integer, intent(in) :: num
real(8), dimension(2), intent(out) :: pt

if(g_cell%g_type.ne.'xy ') call error_message

select case(g_cell%edge(num))
 case(1); pt(1)=g_cell%mdis_posi(num); pt(2)=-0.5d0
 case(2); pt(1)=0.5d0; pt(2)=g_cell%mdis_posi(num)
 case(3); pt(1)=g_cell%mdis_posi(num); pt(2)=0.5d0
 case(4); pt(1)=-0.5d0; pt(2)=g_cell%mdis_posi(num)
end select

end subroutine pick_middle_point


subroutine pick_edge(num,pt1,pt2)
! This subroutine picks the two points of the 'num'th cell-edge
! of a grid cell.

implicit none
real(8), dimension(2), intent(out) :: pt1, pt2
integer, intent(in) :: num

select case(num)
 case(1)
  pt1(1)=-0.5d0; pt1(2)=-0.5d0; pt2(1)=0.5d0; pt2(2)=-0.5d0 
 case(2)
  pt1(1)=0.5d0; pt1(2)=-0.5d0; pt2(1)=0.5d0; pt2(2)=0.5d0
 case(3)
  pt1(1)=0.5d0; pt1(2)=0.5d0; pt2(1)=-0.5d0; pt2(2)=0.5d0
 case(4)
  pt1(1)=-0.5d0; pt1(2)=0.5d0; pt2(1)=-0.5d0; pt2(2)=-0.5d0
end select

end subroutine pick_edge


subroutine find_single(g_cell,single)

! This subroutine finds the single point for 'xy' type critical 
! cell that is on a side of the discontinuity.

implicit none
type(geo_info), intent(in) :: g_cell
integer, intent(out) :: single

if(g_cell%g_type.ne.'xy ') then
 print*, 'Something is wrong!!!'; pause
end if

if(g_cell%point(1).eq.g_cell%point(2)) then
 if(g_cell%point(1).eq.g_cell%point(4)) then
  single=3
 else
  single=4
 end if
else
 if(g_cell%point(1).eq.g_cell%point(4)) then
  single=2
 else
  single=1
 end if
end if

end subroutine find_single


subroutine pick_corner(pt,num)
! This subroutine picks the point 'pt' corresponding to the 
! 'num'th corner of a grid cell.

implicit none
real(8), dimension(2), intent(out) :: pt
integer, intent(in) :: num

select case(num)
 case(1); pt(1)=-0.5d0; pt(2)=-0.5d0
 case(2); pt(1)=0.5d0; pt(2)=-0.5d0
 case(3); pt(1)=0.5d0; pt(2)=0.5d0
 case(4); pt(1)=-0.5d0; pt(2)=0.5d0
end select

end subroutine pick_corner


function side_area(g_cell,side)	result(c)
! This subroutine finds the area of an assigned side.

implicit none
type(geo_info), intent(in) :: g_cell
character*6, intent(in) :: side
real(8) ::c

real(8), dimension(:,:), allocatable :: pt
real(8), dimension(5,2) :: pt1
real(8), dimension(2) :: ptt1, ptt2
real(8) :: area_l, area_r, area
integer :: i, vortices, edge, edge_1, edge_2

area_l=error_data; area_r=error_data
pt1=error_data

call pick_point(g_cell,1,ptt1); pt1(1,:)=ptt1
vortices=1
edge_1=g_cell%edge(1); edge_2=g_cell%edge(2)
if(edge_2.lt.edge_1) edge_2=edge_2+4
do i=edge_1+1, edge_2
 edge=cycle_4(i); call pick_edge(edge,ptt1,ptt2)
 vortices=vortices+1; pt1(vortices,:)=ptt1
end do
vortices=vortices+1; call pick_point(g_cell,2,ptt2)
pt1(vortices,:)=ptt2

allocate(pt(vortices,2))
do i=1, vortices
 pt(i,:)=pt1(i,:)
end do
area=area_of_polygon(pt,vortices)
edge_2=cycle_4(edge_2)
select case(g_cell%point(edge_2))
 case('left  '); area_l=area; area_r=1.0d0-area
 case('right '); area_r=area; area_l=1.0d0-area
 case default; call error_message
end select
select case(side)
 case('left  '); c=area_l
 case('right '); c=area_r
 case default; call error_message
end select

deallocate(pt)  

end function side_area


subroutine compute_curve_length(g_cell,lenght)
! This subroutine computes the length of the curve segment in the critical cell.

implicit none
type(geo_info), intent(in) :: g_cell
real(8), intent(out) :: lenght

real(8), dimension(2) :: pt1, pt2

call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
lenght=length_of_segment(pt1,pt2)

end subroutine compute_curve_length


subroutine which_side(g_cell1,g_cell2,side)
! Given two critical cells 'g_cell1' and 'g_cell2', this subroutine
! find out on which side of the first critical cellthe second one 
! is. The possible values of the output 'side' are 'left', right'
! or 'known'.

implicit none
type(geo_info), intent(in) :: g_cell1, g_cell2
character*6, intent(out) :: side

integer :: single1, single2

select case(g_cell1%g_type)
  
 case('xx ')
  select case(g_cell2%g_type) 
   case('xx '); side='unknow'
   case('yy '); call error_message
   case('xy ')
    call find_single(g_cell2,single2)
    select case(single2)
     case(1,4)
      if(g_cell1%edge(1).eq.1) then
       side='left  '
      else
       side='right '
      end if
     case(2,3)
      if(g_cell1%edge(1).eq.1) then
       side='right '
      else
       side='left  '
      end if
    end select
  end select
  
 case('yy ')
  select case(g_cell2%g_type)
   case('xx '); call error_message
   case('yy '); side='unknow'
   case('xy ')
    call find_single(g_cell2,single2)
    select case(single2)
     case(1,2)
      if(g_cell1%edge(1).eq.2) then
       side='left  '
      else
       side='right '
      end if
     case(3,4)
      if(g_cell1%edge(1).eq.4) then
       side='right '
      else
       side='left  '
      end if
    end select
  end select

 case('xy ')
  call find_single(g_cell1,single1)
  select case(single1)
     
   case(1)
    if(g_cell1%edge(1).eq.1) then
     side='right '
    else
     side='left  '
    end if
      
   case(2)
    if(g_cell1%edge(1).eq.1) then
     side='left  '
    else
     side='right '
    end if
     
   case(3)
    if(g_cell1%edge(1).eq.2) then
     side='left  '
    else
     side='right '
    end if
    
   case(4)
    if(g_cell1%edge(1).eq.3) then
     side='left  '
    else
     side='right '
    end if
  end select 
   
  if(g_cell2%g_type.eq.'xy ') then
   call find_single(g_cell2,single2)
   if(single1.eq.single2) side='unknown'
  end if

end select     	  	 	  	 	 

end subroutine which_side


subroutine normal_4_isolated_xy(g_cell,normal)

implicit none

type(geo_info), intent(in) :: g_cell
real(8), dimension(2) :: normal

real(8), dimension(2) :: pt1, pt2
real(8) :: rr
integer :: single

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
rr=length_of_segment(pt1,pt2)

if(rr.lt.1.0d-2) then
 call find_single(g_cell,single)
 select case(single)
  case(1)
   if(g_cell%edge(1).eq.1) then
    pt1=(/1.0d0,0.0d0/); pt2=(/0.0d0,1.0d0/) 
   else
    pt1=(/0.0d0,1.0d0/); pt2=(/1.0d0,0.0d0/)
   end if
  case(2)
   if(g_cell%edge(1).eq.1) then
    pt1=(/-1.0d0,0.0d0/); pt2=(/0.0d0,1.0d0/) 
   else
    pt1=(/0.0d0,1.0d0/); pt2=(/-1.0d0,0.0d0/)
   end if
  case(3)
   if(g_cell%edge(1).eq.3) then
    pt1=(/-1.0d0,0.0d0/); pt2=(/0.0d0,-1.0d0/) 
   else
    pt1=(/0.0d0,-1.0d0/); pt2=(/-1.0d0,0.0d0/)
   end if
  case(4)
   if(g_cell%edge(1).eq.3) then
    pt1=(/1.0d0,0.0d0/); pt2=(/0.0d0,-1.0d0/) 
   else
    pt1=(/0.0d0,-1.0d0/); pt2=(/1.0d0,0.0d0/)
   end if
 end select 	 
end if

normal=normal_of_line(pt1,pt2)

end subroutine normal_4_isolated_xy


function side_of_point(g_cell,point) result(c)
! This subroutine locates the side the given point is in.

implicit none

type(geo_info), intent(in) :: g_cell
real(8), dimension(2), intent(in) :: point
character*6 :: c

real(8), dimension(2) :: pt1, pt2, normal
real(8) :: distance

call pick_point(g_cell,1, pt1)
call pick_point(g_cell,2, pt2)

normal=normal_of_line(pt1,pt2)
distance=normal(1)*(point(1)-pt1(1))+normal(2)*(point(2)-pt1(2))
if(distance.lt.0.0d0) then
 c='left  '
else
 c='right '
end if

end function side_of_point


function length_of_curve_segment(g_cell) result(c)

implicit none

type(geo_info), intent(in) :: g_cell
real(8) :: c

real(8), dimension(2) :: pt1, pt2

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)

c=length_of_segment(pt1,pt2)

end function length_of_curve_segment


end module geo_info_cell