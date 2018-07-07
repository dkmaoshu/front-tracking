module wind_directions
! This module is for computing the wind directions of the 
! discontinuity curve on the two cell-edges it crosses.

use	auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell, gp_cell, gn_cell
type(phy_info) :: p_cell, pp_cell, pn_cell

private temp, g_cell, gn_cell, p_cell, pn_cell, check_winds
public compute_winds


contains


subroutine compute_winds
! This subroutine computes wind directions on auxiliary 
! discontinuity curves.

implicit none

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call comp_wnds(acvv(i))
  call check_winds(acvv(i))
 end if
end do

end subroutine compute_winds


subroutine comp_wnds(acv)
! Computation on an auxiliary curve.

implicit none

type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, cv_nmb
logical :: head_switch
real(8), dimension(2) :: pt1, pt2, normal, ppt
type(state) :: a, b, c1, c2
real(8) :: speed, velo_x, velo_y
character*12 :: wave_type
integer :: i
real(8) :: x_posi, y_posi

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 address=temp%address; cv_nmb=address%cv_nmb
 g_cell=temp%g_cell; p_cell=temp%p_cell
 do i=1,2
! First, compute the normal at the discontinuity position.

  call clean_up_g_cell(gp_cell); call clean_up_g_cell(gn_cell)
  call clean_up_p_cell(pp_cell); call clean_up_p_cell(pn_cell)
  pt1=error_data; pt2=error_data; ppt=error_data
  x_posi=error_data; y_posi=error_data

  select case(i)

   case(1)
    if(temp%partner.eq.'previous') then
     gn_cell=temp%previous%g_cell; pn_cell=temp%previous%p_cell
    else
     gn_cell=temp%g_cell; pn_cell=temp%p_cell
    end if
    if(associated(temp%p_nxt_dr)) then
     if(temp%p_nxt_dr%partner.eq.'next    ') then
      gp_cell=temp%p_nxt_dr%next%g_cell
      pp_cell=temp%p_nxt_dr%next%p_cell
     else
      gp_cell=temp%p_nxt_dr%g_cell; pp_cell=temp%p_nxt_dr%p_cell
     end if
    else
     gp_cell=gn_cell; pp_cell=pn_cell
    end if
	call pick_point(g_cell,1,ppt)

   case(2)
    if(temp%partner.eq.'next    ') then
     gp_cell=temp%next%g_cell; pp_cell=temp%next%p_cell
    else
     gp_cell=temp%g_cell; pp_cell=temp%p_cell
    end if
    if(associated(temp%n_nxt_dr)) then
     if(temp%n_nxt_dr%partner.eq.'previous') then
      gn_cell=temp%n_nxt_dr%previous%g_cell
      pn_cell=temp%n_nxt_dr%previous%p_cell
     else
      gn_cell=temp%n_nxt_dr%g_cell; pn_cell=temp%n_nxt_dr%p_cell
     end if
    else
     gn_cell=gp_cell; pn_cell=pp_cell
    end if
	call pick_point(g_cell,2,ppt)

  end select

  call pick_point(gp_cell,1,pt1); call pick_point(gn_cell,2,pt2)
  pt2=pt2+(/dfloat(gn_cell%x_idx-gp_cell%x_idx), &
            dfloat(gn_cell%y_idx-gp_cell%y_idx)/)
  normal=normal_of_line(pt1,pt2)
  x_posi=(ppt(1)+dfloat(g_cell%x_idx))*h
  y_posi=(ppt(2)+dfloat(g_cell%y_idx))*h
  
! Then, compute the normal velocity of the discontinuity curve.
  a=0.5d0*(pp_cell%l_state+pn_cell%l_state)
  b=0.5d0*(pp_cell%r_state+pn_cell%r_state)
  call riemann(a,b,normal,p_cell%wv_nb,x_posi,y_posi,c1,c2,speed,wave_type)
  velo_x=dir_velo(c1,c2,f); velo_y=dir_velo(c1,c2,g)

! Finally, determine the wind direction according to the normal 
! velocity.
  select case(g_cell%edge(i))
   case(1)
    if(velo_y.gt.0.d0) then
     g_cell%wind(i)='in '
    else
     g_cell%wind(i)='out'
    end if
   case(2)
    if(velo_x.gt.0.d0) then
     g_cell%wind(i)='out'
    else
     g_cell%wind(i)='in'
    end if
   case(3)
    if(velo_y.gt.0.d0) then
     g_cell%wind(i)='out'
    else
     g_cell%wind(i)='in'
    end if
   case(4)
    if(velo_x.gt.0.d0) then
     g_cell%wind(i)='in '
    else
     g_cell%wind(i)='out'
    end if
  end select

 end do

 if(.not.associated(temp%p_nxt_dr)) then
  if(acv%begin_end.le.0) g_cell%wind(1)='out'
 end if
 if(.not.associated(temp%n_nxt_dr)) then
  if(acv%end_end.le.0) g_cell%wind(2)='out'
 end if

 temp%g_cell=g_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine comp_wnds


subroutine check_winds(acv)
! This subroutine is written for checking whether the winds in
! neighboring critical cells are consistent.

implicit none
type(auxi_discv), intent(inout) :: acv

integer :: head_mark
logical :: head_switch

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch)
 g_cell=temp%g_cell

! Check the consistency of the first wind.
 if(g_cell%wind(1).ne.'in '.and.g_cell%wind(1).ne.'out') then
  call error_message
 end if
 if(associated(temp%p_nxt_dr)) then
  gn_cell=temp%p_nxt_dr%g_cell
  if(gn_cell%wind(2).ne.'in '.and.gn_cell%wind(2).ne.'out') then
   call error_message
  end if
  if(g_cell%wind(1).eq.gn_cell%wind(2)) call error_message
 end if
  
! Check the consistency of the second wind.
 if(g_cell%wind(2).ne.'in '.and.g_cell%wind(2).ne.'out') then
  call error_message
 end if
 if(associated(temp%n_nxt_dr)) then
  gn_cell=temp%n_nxt_dr%g_cell
  if(gn_cell%wind(1).ne.'in '.and.gn_cell%wind(1).ne.'out') then
   call error_message
  end if
  if(g_cell%wind(2).eq.gn_cell%wind(1)) call error_message
 end if

! Check the consistency of winds with the partner.
 if(temp%partner.ne.'single  ') then
  select case(temp%partner)
   case('previous')
    gn_cell=temp%previous%g_cell
   case('next    ')
    gn_cell=temp%next%g_cell
  end select
  if(gn_cell%wind(1).ne.g_cell%wind(1).or.gn_cell% &
     wind(2).ne.g_cell%wind(2)) call error_message
 end if  	   

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine check_winds


end module wind_directions