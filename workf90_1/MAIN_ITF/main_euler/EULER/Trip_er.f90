module tools_used_in_triple

use physical_state
! 'euler.f90'

use euler_functions
! 'euler_functions.f90'

use tools
!' tools.f90'

implicit none

character*6, dimension(3) :: wave_facing=(/'left  ', 'no    ', 'right '/)
! To indicating whether the discontinuity wave is facing to the
! left, right or not specified. 

public  triple_multiple, decompose_conservation


contains


subroutine triple_multiple(curves_map,tri_multiple)
! The final physical check whether a node cell is triple or 
! multiple is done by this subroutine based on the information
! provided by the subroutine 'check_triple' in module 
! 'connect_and_triple' coded in Fortran90 file 'contrip.f90'.
! The performance of this subroutine is really physically related;
! therefore, we put it in this physical subdirectory to make the main
! body of the code, precisely speaking the module 'connect_and_ 
! triple', be physically independent, which is necessary for
! maintaining the code and extending it to the Euler system later.

implicit none

type(curve_map), dimension(wave_number), intent(inout) ::  &
curves_map
! Tells wheter the discontinuity is facing to the previous or 
! next smooth region and its status.
character*8, intent(out) :: tri_multiple
! To tell whether the node cell is triple or multiple.

! The configuration of triple points for Euler equations is as follows: 
! The 'oncoming' discontinuities are at most of two, the two 'outgoing' shocks are
! on the two sides with the left one going to the left and the right one the right,
! the contact discontinuity is in between. One of the 'outgoing' discontinuities may
! be missing, but never two.

integer :: i, out_number

! Compute the number of outgoing discontinuities.
!out_number=out_number+1
out_number=0
do i=1,wave_number
 if(curves_map(i)%status.eq.'awake '.or.curves_map(i)%status.eq.'yawn  ') then
  out_number=out_number+1
 end if
end do

select case(out_number)

 case(0); tri_multiple='remove  '
 
 case(1); tri_multiple='triple  '
!  if(curves_map(1)%status.eq.'yawn  ') then
!   tri_multiple='triple  '
!  else
!   tri_multiple='triple  '
!  end if
!  return

 case(2)
  if(curves_map(1)%status.eq.'yawn  '.or.curves_map(2)%status.eq.'yawn  ') then
   tri_multiple='mulitiple'; return
  end if
  if(curves_map(1)%facing.eq.'next    ') then
   tri_multiple='multiple'; return
  end if
  if(curves_map(1)%facing.eq.'previous') then
   if(curves_map(2)%facing.ne.'previous') then
    tri_multiple='triple  '; return
   else
    tri_multiple='multiple'; return
   end if
  end if
  if(curves_map(1)%facing.eq.'no      ') then
   if(curves_map(2)%facing.eq.'next    ') then
    tri_multiple='triple  '; return
   else
    tri_multiple='multiple'; return
   end if
  end if

 case(3)
  do i=1,3
   if(curves_map(i)%status.eq.'yawn  ') then
    tri_multiple='multiple'
   end if
  end do
  if(curves_map(i)%facing.ne.'previous') then
   tri_multiple='multiple'; return
  end if
  if(curves_map(2)%facing.ne.'no      ') then
   tri_multiple='multiple'; return
  end if
  tri_multiple='triple  '; return
   	     	 
end select

end subroutine triple_multiple


subroutine decompose_conservation(l_state,r_state,normal,sum,wave_number,dist_conv)
! This subroutine decomposes a given quantity based on Riemann
! separation between left and right states on a normal direction.

implicit none
type(state), intent(in) :: l_state, r_state
real(8), dimension(2), intent(in) :: normal
integer, intent(in) :: wave_number
type(state), intent(in) :: sum
type(state), dimension(:), intent(out) :: dist_conv

type(state) :: l_state_in_normal, r_state_in_normal
type(state) :: c1, c2, c1_in_normal, c2_in_normal
type(state) :: intermediate, intermediate_in_normal
type(state), dimension(4) :: left_vectors_in_normal, right_vectors_in_normal
type(state) :: right_vector
real(8) :: velocity, check
character*12 :: wave_type
real(8), dimension(4) :: x, sma, bb

real(8), dimension(4,4) :: qq
integer :: size_of_dc, j

real(8) :: x_temp, y_temp

size_of_dc=size(dist_conv)
if(size_of_dc.ne.3) print*, 'There must be something wrong!'

x_temp=0.0d0; y_temp=0.0d0
call riemann(l_state,r_state,normal,2,x_temp,y_temp,c1,c2,velocity, &
             wave_type)

select case(wave_number)

 case(1)
  do j=1, 4
   qq(j,1)=c1%value(j)-l_state%value(j)
  end do
  call transform_1(c1,normal,c1_in_normal)
  call transform_1(c2,normal,c2_in_normal)
  call roe_decomposition(c1_in_normal,c2_in_normal,left_vectors_in_normal, &
                         right_vectors_in_normal)
  call transform_2(right_vectors_in_normal(2),normal,right_vector)
  do j=1, 4
   qq(j,2)=right_vector%value(j)
  end do
  call transform_2(right_vectors_in_normal(3),normal,right_vector)
  do j=1, 4
   qq(j,3)=right_vector%value(j)
  end do
  call transform_1(r_state,normal,r_state_in_normal)
  call roe_decomposition(c2_in_normal,r_state_in_normal,left_vectors_in_normal, &
                         right_vectors_in_normal)
  call transform_2(right_vectors_in_normal(4),normal,right_vector)
  do j=1, 4
   qq(j,4)=right_vector%value(j)
  end do
 

 case(3)
  do j=1, 4
   qq(j,4)=r_state%value(j)-c2%value(j)
  end do
  call transform_1(c1,normal,c1_in_normal)
  call transform_1(c2,normal,c2_in_normal)
  call roe_decomposition(c1_in_normal,c2_in_normal,left_vectors_in_normal, &
                         right_vectors_in_normal)
  call transform_2(right_vectors_in_normal(2),normal,right_vector)
  do j=1, 4
   qq(j,2)=right_vector%value(j)
  end do
  call transform_2(right_vectors_in_normal(3),normal,right_vector)
  do j=1, 4
   qq(j,3)=right_vector%value(j)
  end do
  call transform_1(l_state,normal,l_state_in_normal)
  call roe_decomposition(l_state_in_normal,c1_in_normal,left_vectors_in_normal, &
                         right_vectors_in_normal)
  call transform_2(right_vectors_in_normal(1),normal,right_vector)
  do j=1, 4
   qq(j,1)=right_vector%value(j)
  end do

 case(2)
  call transform_1(l_state,normal,l_state_in_normal)
  call transform_1(c1,normal,c1_in_normal)
  call roe_decomposition(l_state_in_normal,c1_in_normal,left_vectors_in_normal, &
                         right_vectors_in_normal)
  call transform_2(right_vectors_in_normal(1),normal,right_vector)
  do j=1, 4
   qq(j,1)=right_vector%value(j)
  end do
  call transform_1(c2,normal,c2_in_normal)
  call transform_1(r_state,normal,r_state_in_normal)
  call roe_decomposition(c2_in_normal,r_state_in_normal,left_vectors_in_normal, &
                         right_vectors_in_normal)
  call transform_2(right_vectors_in_normal(4), normal,right_vector)
  do j=1, 4
   qq(j,4)=right_vector%value(j)
  end do
  call form_intermediate(c1_in_normal,c2_in_normal,intermediate_in_normal) 
  call transform_2(intermediate_in_normal,normal,intermediate)
  do j=1, 4
   qq(j,2)=intermediate%value(j)-c1%value(j)
   qq(j,3)=c2%value(j)-intermediate%value(j)
  end do
  check=dmax1(dabs(qq(2,1)),dabs(qq(2,2)),dabs(qq(2,3)),dabs(qq(2,4)))
  if(check.lt.1.0d-3) then
   call transform_1(c1,normal,c1_in_normal)
   call transform_1(c2,normal,c2_in_normal)
   call roe_decomposition(c1_in_normal,c2_in_normal,left_vectors_in_normal, &
                          right_vectors_in_normal)
   call transform_2(right_vectors_in_normal(2),normal,right_vector)
   do j=1, 4
    qq(j,2)=right_vector%value(j)
   end do
  end if
  check=dmax1(dabs(qq(3,1)),dabs(qq(3,2)),dabs(qq(3,3)),dabs(qq(3,4)))
  if(check.lt.1.0d-3) then
   call transform_1(c1,normal,c1_in_normal)
   call transform_1(c2,normal,c2_in_normal)
   call roe_decomposition(c1_in_normal,c2_in_normal,left_vectors_in_normal, &
                          right_vectors_in_normal)
   call transform_2(right_vectors_in_normal(3),normal,right_vector)
   do j=1, 4
    qq(j,3)=right_vector%value(j)
   end do
  end if

end select
call pivots_deletion(qq,sum%value,x,sma)
do j=1, 4
 dist_conv(1)%value(j)=x(4)*qq(j,4)
 dist_conv(2)%value(j)=x(2)*qq(j,2)+x(3)*qq(j,3)
 dist_conv(3)%value(j)=x(1)*qq(j,1)
end do


contains


subroutine form_intermediate(left_state,right_state,intermediate_state)

implicit none

type(state), intent(in) :: left_state, right_state
type(state), intent(out) :: intermediate_state
real(8), dimension(4) :: wwl, wwr, ww
real(8) :: rhol, unl, utl, pl, rhor, unr, utr, pr, rho, un, ut, p  

wwl=left_state%value
rhol=wwl(1)
unl=uf(wwl)
utl=vf(wwl)
pl=pf(wwl,left_state%gamma)
wwr=right_state%value
rhor=wwr(1)
unr=uf(wwr)
utr=vf(wwr)
pr=pf(wwr,right_state%gamma)
rho=rhol
un=unr
ut=utr
p=pr
call tran3(rho,un,ut,p,ww,left_state%gamma)
intermediate_state%value=ww
intermediate_state%gamma=left_state%gamma

end subroutine form_intermediate


subroutine transform_1(given_state,normal,transformed_state)
! This subroutine transforms the 'given_state' in 'xy'-coordinates into 
! the 'transformed_state' in 'normal-tangential'-coordinates.

implicit none
type(state), intent(in) :: given_state
real(8), dimension(2) :: normal
type(state), intent(out) :: transformed_state

real(8), dimension(2) :: mxy, mnt, tangent

tangent(1)=-normal(2); tangent(2)=normal(1)
transformed_state=given_state
mxy(1)=given_state%value(2); mxy(2)=given_state%value(3)
mnt(1)=inner_product(normal,mxy)
mnt(2)=inner_product(tangent,mxy)
transformed_state%value(2)=mnt(1)
transformed_state%value(3)=mnt(2)

end subroutine transform_1


subroutine transform_2(given_state,normal,transformed_state)
! This subroutine transforms the 'given_state' in 'xy'-coordinates into
! the 'transformed_state' in 'normal-tangential' coordinates.

implicit none
type(state), intent(in) :: given_state
real(8), dimension(2) :: normal
type(state), intent(out) :: transformed_state

real(8), dimension(2) :: mxy, mnt, tangent

tangent(1)=-normal(2); tangent(2)=normal(1)
transformed_state=given_state
mnt(1)=given_state%value(2)
mnt(2)=given_state%value(3)
mxy(1)=mnt(1)*normal(1)+mnt(2)*tangent(1)
mxy(2)=mnt(1)*normal(2)+mnt(2)*tangent(2)
transformed_state%value(2)=mxy(1)
transformed_state%value(3)=mxy(2)

end subroutine transform_2


function inner_product(vec1,vec2) result(c)

implicit none
real(8), dimension(:), intent(in) :: vec1, vec2
real(8) :: c

integer :: size1, size2, i

size1=size(vec1); size2=size(vec2)
if(size1.ne.size2) call error_message

c=0.0d0
do i=1,size1
 c=c+vec1(i)*vec2(i)
end do

end function inner_product


subroutine roe_decomposition(state_1,state_2,left_eigenvectors,right_eigenvectors)
! Given physical states 'state_1' and 'state_2', finds the eigenvectors of
! the corresponding Roe's Riemann approximation solver.

implicit none
type(state), intent(in) :: state_1, state_2
type(state), dimension(4), intent(out) :: left_eigenvectors, right_eigenvectors

real(8), dimension(4) :: w1, w2, ch, coef
real(8), dimension(4,4) :: el, er
integer :: i

w1=state_1%value; w2=state_2%value
call roe(w2,w1,el,er,ch,coef,state_1%gamma)
do i=1,4
 left_eigenvectors(i)%value=el(i,:)
 right_eigenvectors(i)%value=er(i,:)
end do

end subroutine roe_decomposition


end subroutine decompose_conservation


end module tools_used_in_triple
