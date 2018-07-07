# Microsoft Developer Studio Project File - Name="main_euler" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=main_euler - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "main_euler.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "main_euler.mak" CFG="main_euler - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "main_euler - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "main_euler - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "main_euler - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir ".\Release"
# PROP BASE Intermediate_Dir ".\Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\Release"
# PROP Intermediate_Dir ".\Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /I "Release/"
# ADD F90 /compile_only /include:"Release/" /math_library:fast /nologo
# ADD CPP /FD
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir ".\Debug"
# PROP BASE Intermediate_Dir ".\Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\Debug"
# PROP Intermediate_Dir ".\Debug"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /debug:full /nologo /I "Debug/"
# ADD F90 /browser /compile_only /debug:full /include:"Debug/" /nologo
# ADD CPP /FR /FD
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386

!ENDIF 

# Begin Target

# Name "main_euler - Win32 Release"
# Name "main_euler - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;for;f90"
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\adss_cell.f90
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\auxi_cuv\auxi_cc.f90
NODEP_F90_AUXI_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\computational_information.mod"\
	".\Debug\critical_cells.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\prc_auxi\auxi_ccp.f90
NODEP_F90_AUXI_C=\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\AUXI_CUV\Auxi_cv1.f90
NODEP_F90_AUXI_CV=\
	".\Debug\auxiliary_critical_cells.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\prc_auxi\auxi_cvp.f90
NODEP_F90_AUXI_CVP=\
	".\Debug\auxi_cc_prod.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\auxiliary_list_fix.mod"\
	".\Debug\clean_up.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\node_cells.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\boundary\bd_posi.f90
NODEP_F90_BD_PO=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\boundary\boundary_conditions.f90
DEP_F90_BOUND=\
	".\Debug\Haas_sturtevant_1d_boundary.mod"\
	".\Debug\Haas_Sturtevant_boundary.mod"\
	".\Debug\Jacobs_Peng_Zabusky_boundary.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\bs_var.f90
NODEP_F90_BS_VA=\
	".\Debug\polygon_intersection.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\chk_rmv.f90
NODEP_F90_CHK_R=\
	".\Debug\tools_4_connecting_n_triple.mod"\
	".\Debug\tools_used_in_triple.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PHYS_FLX\PHY_UPDT\Clean_up6.f90
NODEP_F90_CLEAN=\
	".\Debug\discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\phys_flx\phy_updt\comlids.f90
NODEP_F90_COMLI=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_can.f90
NODEP_F90_COMP_=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\geo_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\comp_circ_HS.f90
NODEP_F90_COMP_C=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_dpn.f90
NODEP_F90_COMP_D=\
	".\Debug\pick_middle_positions.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_sth\comp_in_smooth.f90
NODEP_F90_COMP_I=\
	".\Debug\boundary_conditions.mod"\
	".\Debug\node_cells.mod"\
	".\Debug\numerical_fluxes.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_nd.f90
NODEP_F90_COMP_N=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\grid_map.mod"\
	".\Debug\node_cells.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_on_cvn.f90
NODEP_F90_COMP_O=\
	".\Debug\boundary_conditions.mod"\
	".\Debug\compute_boundary_positions.mod"\
	".\Debug\compute_cell_average.mod"\
	".\Debug\compute_discontinuity_positions.mod"\
	".\Debug\compute_node_ca.mod"\
	".\Debug\half_step_computation.mod"\
	".\Debug\numerical_regularization.mod"\
	".\Debug\output_middle_positions.mod"\
	".\Debug\recover_mid_positions.mod"\
	".\Debug\restriction_of_positions.mod"\
	".\Debug\update_left_n_right_states.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\compute_growth_rate.f90
NODEP_F90_COMPU=\
	".\Debug\discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\contrip.f90
NODEP_F90_CONTR=\
	".\Debug\check_n_remove.mod"\
	".\Debug\fix_broken_connection.mod"\
	".\Debug\remove_broken_corners.mod"\
	".\Debug\treatment_of_triples.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\reg_cuv\crit_cell.f90
NODEP_F90_CRIT_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\crit_info.f90
NODEP_F90_CRIT_I=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\solu\data_ext.f90
NODEP_F90_DATA_=\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\data_pre.f90
NODEP_F90_DATA_P=\
	".\Debug\node_cells.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\interact\dis_jntn.f90
NODEP_F90_DIS_J=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\reg_cuv\discontinuity_curve.f90
NODEP_F90_DISCO=\
	".\Debug\critical_cells.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\euler\euler_functions.f90
# End Source File
# Begin Source File

SOURCE=.\EULER\eulernew.f90
NODEP_F90_EULER=\
	".\Debug\data_extrapolation.mod"\
	".\Debug\euler_functions.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\phys_flx\phy_updt\eva_sts.f90
NODEP_F90_EVA_S=\
	".\Debug\euler_functions.mod"\
	".\Debug\produce_polynomial_4_xx_yy.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\auxi_cuv\fif_cell.f90
NODEP_F90_FIF_C=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\fix_cnnt.f90
NODEP_F90_FIX_C=\
	".\Debug\tools_4_connecting_n_triple.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\fix_recv.f90
NODEP_F90_FIX_R=\
	".\Debug\node_cells.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\fix_trip.f90
NODEP_F90_FIX_T=\
	".\Debug\tools_4_connecting_n_triple.mod"\
	".\Debug\tools_used_in_triple.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PHYS_FLX\weno\flux_function_weno2.f90
NODEP_F90_FLUX_=\
	".\Debug\euler_functions.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PHYS_FLX\weno\flux_weno.f90
NODEP_F90_FLUX_W=\
	".\Debug\flux_functions.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\form_dsv.f90
NODEP_F90_FORM_=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\fun_use.f90
NODEP_F90_FUN_U=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\gif_cell.f90
NODEP_F90_GIF_C=\
	".\Debug\grid.mod"\
	".\Debug\tools.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\grid.f90
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\grid_map.f90
NODEP_F90_GRID_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\grid.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\H_stp_eul.f90
NODEP_F90_H_STP=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\INI_DATA\init_HS_airHe\HS1d_bond_cond.f90
NODEP_F90_HS1D_=\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\INI_DATA\init_HS_airHe\HS_boun_cond.f90
NODEP_F90_HS_BO=\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\inputs_euler1.f90
NODEP_F90_INPUT=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\tools\interpl.f90
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\prc_auxi\isl_fix.f90
NODEP_F90_ISL_F=\
	".\Debug\auxi_cc_prod.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\INI_DATA\init_RM_peng\Jacobs_boun_cond1.f90
NODEP_F90_JACOB=\
	".\Debug\numerical_fluxes.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\locl_fix.f90
NODEP_F90_LOCL_=\
	".\Debug\solu_comput.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\main.f90
NODEP_F90_MAIN_=\
	".\Debug\auxiliary_curves_production.mod"\
	".\Debug\compress_linear_dis.mod"\
	".\Debug\computation_in_smooth.mod"\
	".\Debug\computation_on_curves.mod"\
	".\Debug\compute_circulation_HS.mod"\
	".\Debug\input_setup.mod"\
	".\Debug\mesh_ratio.mod"\
	".\Debug\move_all.mod"\
	".\Debug\node_move_n_split.mod"\
	".\Debug\out_intermediate.mod"\
	".\Debug\output_show.mod"\
	".\Debug\output_show_HR.mod"\
	".\Debug\reset_discvs.mod"\
	".\Debug\revers_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\rest_cv\manipulate\manipu_3.f90
NODEP_F90_MANIP=\
	".\Debug\discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\rest_cv\manipulate\manipu_4.f90
NODEP_F90_MANIPU=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\merge_ck.f90
NODEP_F90_MERGE=\
	".\Debug\variables_4_node_treat.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\mesh_ratio\mesh_rt.f90
NODEP_F90_MESH_=\
	".\Debug\numerical_fluxes.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\boundary\mv_all.f90
NODEP_F90_MV_AL=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\node\nd_cls.f90
NODEP_F90_ND_CL=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\node_plug_ring.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\nd_mvsp.f90
NODEP_F90_ND_MV=\
	".\Debug\connect_and_triple.mod"\
	".\Debug\corner_remove.mod"\
	".\Debug\disjoints_fix.mod"\
	".\Debug\merge_checks.mod"\
	".\Debug\node_production_n_update_1.mod"\
	".\Debug\node_production_n_update_2.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_proc1.f90
NODEP_F90_ND_PR=\
	".\Debug\variables_4_node_treat.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_proc2.f90
NODEP_F90_ND_PRO=\
	".\Debug\variables_4_node_treat.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\node\nd_ring.f90
NODEP_F90_ND_RI=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PHYS_FLX\PHY_UPDT\NS_diffusion.f90
NODEP_F90_NS_DI=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\num_int.f90
NODEP_F90_NUM_I=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\num_int2.f90
NODEP_F90_NUM_IN=\
	".\Debug\function_in_use.mod"\
	".\Debug\numerical_integrals.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\COMP_CV\num_regular5.f90
NODEP_F90_NUM_R=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\geo_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\REST_CV\num_regular_partn5.f90
NODEP_F90_NUM_RE=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\geo_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\out_mdl_euler.f90
NODEP_F90_OUT_M=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\output_for_show_euler.f90
NODEP_F90_OUTPU=\
	".\Debug\physical_state.mod"\
	".\Debug\scan_conservation_on_list.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\outputs_euler1.f90
NODEP_F90_OUTPUT=\
	".\Debug\grid_map.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\REST_CV\packed_separate.f90
NODEP_F90_PACKE=\
	".\Debug\discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\pick_dp.f90
NODEP_F90_PICK_=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\pick_mdp.f90
NODEP_F90_PICK_M=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\pif_cell.f90
NODEP_F90_PIF_C=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\ploynml.f90
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\polyin.f90
NODEP_F90_POLYI=\
	".\Debug\tools.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\interact\pro_lcl.f90
NODEP_F90_PRO_L=\
	".\Debug\solu_comput.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\pro_smnd.f90
NODEP_F90_PRO_S=\
	".\Debug\produce_local_smooth.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\phys_flx\proc_pln.f90
NODEP_F90_PROC_=\
	".\Debug\interpolation.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\COMP_CV\Rcv_mpne3.f90
NODEP_F90_RCV_M=\
	".\Debug\pick_discontinuity_positions.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\reset.f90
NODEP_F90_RESET=\
	".\Debug\boundary_conditions.mod"\
	".\Debug\clean_up.mod"\
	".\Debug\data_preparation.mod"\
	".\Debug\disjoints_fix.mod"\
	".\Debug\fix_reversed_curves.mod"\
	".\Debug\form_discurves.mod"\
	".\Debug\local_smooth_fix.mod"\
	".\Debug\manipulation_4.mod"\
	".\Debug\numerical_NS_diffusion.mod"\
	".\Debug\numerical_regular_partners.mod"\
	".\Debug\packed_separation.mod"\
	".\Debug\produce_local_smooth_4_nodes.mod"\
	".\Debug\smoothen_curves.mod"\
	".\Debug\update_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\restr_dp.f90
NODEP_F90_RESTR=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\rev_cvs.f90
NODEP_F90_REV_C=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\rm_b_crn.f90
NODEP_F90_RM_B_=\
	".\Debug\tools_4_connecting_n_triple.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\rm_cons.f90
NODEP_F90_RM_CO=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\rest_cv\rmv_tgls.f90
NODEP_F90_RMV_T=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\solu_comput.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\save_growth_rate.f90
NODEP_F90_SAVE_=\
	".\Debug\compute_growth_rate.mod"\
	".\Debug\input_setup.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\scan\sc_conls.f90
NODEP_F90_SC_CO=\
	".\Debug\discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\euler\sh_sths.f90
NODEP_F90_SH_ST=\
	".\Debug\euler_functions.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\interface\sh_sts_el_HR.f90
NODEP_F90_SH_STS=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\REST_CV\smooth_curves3.f90
NODEP_F90_SMOOT=\
	".\Debug\discontinuity_curves.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\solu\solu_com.f90
NODEP_F90_SOLU_=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\grid_map.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\solu\solution.f90
NODEP_F90_SOLUT=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\grid_map.mod"\
	".\Debug\node_cells.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_strct\sl_smth\solution_smooth.f90
NODEP_F90_SOLUTI=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	".\Debug\show_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\tools.f90
NODEP_F90_TOOLS=\
	".\Debug\polynomials.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\tools_ct.f90
NODEP_F90_TOOLS_=\
	".\Debug\variables_4_node_treat.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\REST_CV\tools_rtn2.f90
NODEP_F90_TOOLS_R=\
	".\Debug\solu_comput.mod"\
	".\Debug\tools.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\euler\trip_er.f90
NODEP_F90_TRIP_=\
	".\Debug\euler_functions.mod"\
	".\Debug\physical_state.mod"\
	".\Debug\tools.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\interact\nd_mvspl\trt_trpl.f90
NODEP_F90_TRT_T=\
	".\Debug\fix_triple_nodes.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\upd_cvs.f90
NODEP_F90_UPD_C=\
	".\Debug\remove_tangles.mod"\
	".\Debug\upd_cv_cases_1.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\phys_flx\phy_updt\updat_lr.f90
NODEP_F90_UPDAT=\
	".\Debug\evaluate_updated_xxyy.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\updcvcas.f90
NODEP_F90_UPDCV=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\numerical_integral_in_2D.mod"\
	".\Debug\solu_comput.mod"\
	".\Debug\tools_for_reset.mod"\
	
# End Source File
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\updcvcas_1n.f90
NODEP_F90_UPDCVC=\
	".\Debug\node_cells.mod"\
	".\Debug\upd_cv_cases.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;cnt;rtf;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
