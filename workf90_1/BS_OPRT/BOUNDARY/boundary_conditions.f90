module boundary_conditions

!use richtmyer_meshkov_boundary
!! 'rm_bound_cond.f90'

!use free_stream_boundary
!! 'fs_bound_cond.f90'

use Haas_sturtevant_boundary
! 'HS_bound_cond.f90'

use Jacobs_Peng_Zabusky_boundary
! 'Jacobs_boun_cond.f90'

use Haas_Sturtevant_1d_boundary
! 'HS1d_boun_cond.f90'

implicit none


public  boundary_condition1, boundary_condition2, boundary_condition3
private Haas_sturtevant_1, &
		Haas_Sturtevant_2, Haas_Sturtevant_3, Jacobs_PengZ_1, Jacobs_PengZ_2, &
		Jacobs_PengZ_3, Haas_Sturtevant_1d_1, Haas_Sturtevant_1d_2, &
        Haas_Sturtevant_1d_3


contains


subroutine boundary_condition1

implicit none

select case(boundary_type_1)

! case('RM_instabili'); call richtmyer_meshkov_1
! case('free_stream '); call free_stream_1
 case('Haas_Sturtev'); call Haas_sturtevant_1
 case('Jacobs_PengZ'); call Jacobs_PengZ_1
 case('Haas_Stur_1d'); call Haas_Sturtevant_1d_1
 case default; call error_message

end select

end subroutine boundary_condition1


subroutine boundary_condition2

implicit none

select case(boundary_type_1)

! case('RM_instabili'); call richtmyer_meshkov_2
! case('free_stream '); call free_stream_2
 case('Haas_Sturtev'); call Haas_sturtevant_2
 case('Jacobs_PengZ'); call Jacobs_PengZ_2
 case('Haas_Stur_1d'); call Haas_Sturtevant_1d_2
 case default; call error_message

end select

end subroutine boundary_condition2


subroutine boundary_condition3

implicit none

select case(boundary_type_1)

! case('RM_instabili'); call richtmyer_meshkov_3
! case('free_stream '); call free_stream_3
 case('Haas_Sturtev'); call Haas_Sturtevant_3
 case('Jacobs_PengZ'); call Jacobs_PengZ_3
 case('Haas_Stur_1d'); call Haas_Sturtevant_1d_3
 case default; call error_message

end select

end subroutine boundary_condition3


subroutine nonreflecting_boundary_1(uum)

implicit none

type(state), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2), intent(inout) :: uum

select case(boundary_type_1)

 case('Jacobs_PengZ'); call Jacobs_PZ_nonreflecting_1(uum)
 case default; call error_message

end select

end subroutine nonreflecting_boundary_1


subroutine nonreflecting_boundary_2(uum,uu1)

implicit none

type(state), dimension(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2), intent(in) :: uum, uu1

select case(boundary_type_1)

 case('Jacobs_PengZ'); call Jacobs_PZ_nonreflecting_2(uum,uu1)
 case default; call error_message

end select

end subroutine nonreflecting_boundary_2


end module boundary_conditions