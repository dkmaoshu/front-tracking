# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=main_euler - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to main_euler - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "main_euler - Win32 Release" && "$(CFG)" !=\
 "main_euler - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "main_euler.mak" CFG="main_euler - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "main_euler - Win32 Release" (based on\
 "Win32 (x86) Console Application")
!MESSAGE "main_euler - Win32 Debug" (based on\
 "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "main_euler - Win32 Debug"
RSC=rc.exe
F90=fl32.exe

!IF  "$(CFG)" == "main_euler - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\main_euler.exe" "$(OUTDIR)\critical_cell_information.mod"

CLEAN : 
	-@erase ".\Release\critical_cell_information.mod"
	-@erase ".\Release\adss_info_cell.mod"
	-@erase ".\Release\geo_info_cell.mod"
	-@erase ".\Release\grid.mod"
	-@erase ".\Release\tools.mod"
	-@erase ".\Release\polynomials.mod"
	-@erase ".\Release\phy_info_cell.mod"
	-@erase ".\Release\physical_state.mod"
	-@erase ".\Release\euler_functions.mod"
	-@erase ".\Release\data_extrapolation.mod"
	-@erase ".\Release\main_euler.exe"
	-@erase ".\Release\main_temp.obj"
	-@erase ".\Release\input_setup.mod"
	-@erase ".\Release\solution.mod"
	-@erase ".\Release\grid_map.mod"
	-@erase ".\Release\discontinuity_curves.mod"
	-@erase ".\Release\critical_cells.mod"
	-@erase ".\Release\auxiliary_discontinuity_curves.mod"
	-@erase ".\Release\auxiliary_critical_cells.mod"
	-@erase ".\Release\computational_information.mod"
	-@erase ".\Release\solu_in_smth.mod"
	-@erase ".\Release\show_smth.mod"
	-@erase ".\Release\node_cells.mod"
	-@erase ".\Release\node_plug_ring.mod"
	-@erase ".\Release\out_intermediate.mod"
	-@erase ".\Release\output_show.mod"
	-@erase ".\Release\scan_conservation_on_list.mod"
	-@erase ".\Release\auxiliary_curves_production.mod"
	-@erase ".\Release\auxi_cc_prod.mod"
	-@erase ".\Release\smoothen_xy_cell_averages.mod"
	-@erase ".\Release\auxiliary_list_fix.mod"
	-@erase ".\Release\computation_in_smooth.mod"
	-@erase ".\Release\numerical_fluxes.mod"
	-@erase ".\Release\flux_functions.mod"
	-@erase ".\Release\solu_comput.mod"
	-@erase ".\Release\computation_on_curves.mod"
	-@erase ".\Release\compute_cell_average.mod"
	-@erase ".\Release\recover_mid_positions.mod"
	-@erase ".\Release\pick_discontinuity_positions.mod"
	-@erase ".\Release\compute_discontinuity_positions.mod"
	-@erase ".\Release\pick_middle_positions.mod"
	-@erase ".\Release\restriction_of_positions.mod"
	-@erase ".\Release\half_step_computation.mod"
	-@erase ".\Release\compute_node_ca.mod"
	-@erase ".\Release\output_middle_positions.mod"
	-@erase ".\Release\compute_boundary_positions.mod"
	-@erase ".\Release\reset_discvs.mod"
	-@erase ".\Release\update_left_n_right_states.mod"
	-@erase ".\Release\evaluate_updated_xxyy.mod"
	-@erase ".\Release\produce_polynomial_4_xx_yy.mod"
	-@erase ".\Release\interpolation.mod"
	-@erase ".\Release\produce_local_smooth_4_nodes.mod"
	-@erase ".\Release\produce_local_smooth.mod"
	-@erase ".\Release\local_smooth_fix.mod"
	-@erase ".\Release\data_preparation.mod"
	-@erase ".\Release\update_curves.mod"
	-@erase ".\Release\upd_cv_cases_1.mod"
	-@erase ".\Release\upd_cv_cases.mod"
	-@erase ".\Release\numerical_integral_in_2D.mod"
	-@erase ".\Release\numerical_integrals.mod"
	-@erase ".\Release\function_in_use.mod"
	-@erase ".\Release\tools_for_reset.mod"
	-@erase ".\Release\remove_tangles.mod"
	-@erase ".\Release\fix_reversed_curves.mod"
	-@erase ".\Release\form_discurves.mod"
	-@erase ".\Release\disjoints_fix.mod"
	-@erase ".\Release\mesh_ratio.mod"
	-@erase ".\Release\output_show_HR.mod"
	-@erase ".\Release\sh_sths.obj"
	-@erase ".\Release\auxi_cv.obj"
	-@erase ".\Release\comp_can.obj"
	-@erase ".\Release\restr_dp.obj"
	-@erase ".\Release\fix_recv.obj"
	-@erase ".\Release\sc_conls.obj"
	-@erase ".\Release\flux.obj"
	-@erase ".\Release\solu_com.obj"
	-@erase ".\Release\comp_dpn.obj"
	-@erase ".\Release\tools.obj"
	-@erase ".\Release\data_ext.obj"
	-@erase ".\Release\nd_ring.obj"
	-@erase ".\Release\isl_fix.obj"
	-@erase ".\Release\updcvcas_1n.obj"
	-@erase ".\Release\inputs_euler.obj"
	-@erase ".\Release\auxi_cc.obj"
	-@erase ".\Release\solution.obj"
	-@erase ".\Release\proc_pln.obj"
	-@erase ".\Release\comp_nd.obj"
	-@erase ".\Release\tools_rtn.obj"
	-@erase ".\Release\rmv-tgls.obj"
	-@erase ".\Release\pro_lcl.obj"
	-@erase ".\Release\auxi_cvp.obj"
	-@erase ".\Release\adss_cell.obj"
	-@erase ".\Release\reset.obj"
	-@erase ".\Release\grid_map.obj"
	-@erase ".\Release\solution_smooth.obj"
	-@erase ".\Release\flux_functions.obj"
	-@erase ".\Release\euler_functions.obj"
	-@erase ".\Release\smooth_xy.obj"
	-@erase ".\Release\pif_cell.obj"
	-@erase ".\Release\interpl.obj"
	-@erase ".\Release\fun_use.obj"
	-@erase ".\Release\output_for_show_euler.obj"
	-@erase ".\Release\upd_cvs.obj"
	-@erase ".\Release\pick_dp.obj"
	-@erase ".\Release\euler.obj"
	-@erase ".\Release\comp_on_cvn.obj"
	-@erase ".\Release\auxi_ccp.obj"
	-@erase ".\Release\updat_lr.obj"
	-@erase ".\Release\outputs_euler.obj"
	-@erase ".\Release\discontinuity_curve.obj"
	-@erase ".\Release\pro_smnd.obj"
	-@erase ".\Release\comp_in_smooth.obj"
	-@erase ".\Release\half_stp.obj"
	-@erase ".\Release\mesh_rt.obj"
	-@erase ".\Release\eva_sts.obj"
	-@erase ".\Release\updcvcas.obj"
	-@erase ".\Release\pick_mdp.obj"
	-@erase ".\Release\num_int.obj"
	-@erase ".\Release\data_pre.obj"
	-@erase ".\Release\dis_jntn.obj"
	-@erase ".\Release\crit_cell.obj"
	-@erase ".\Release\crit_info.obj"
	-@erase ".\Release\form_dsv.obj"
	-@erase ".\Release\out_mdl_euler.obj"
	-@erase ".\Release\ploynml.obj"
	-@erase ".\Release\num_int2.obj"
	-@erase ".\Release\nd_cls.obj"
	-@erase ".\Release\rcv_mpn.obj"
	-@erase ".\Release\locl_fix.obj"
	-@erase ".\Release\grid.obj"
	-@erase ".\Release\gif_cell.obj"
	-@erase ".\Release\bd_posi.obj"
	-@erase ".\Release\fif_cell.obj"
	-@erase ".\Release\sh_sts_el_HR.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/main_euler.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/main_euler.pdb" /machine:I386 /out:"$(OUTDIR)/main_euler.exe" 
LINK32_OBJS= \
	".\Release\main_temp.obj" \
	".\Release\sh_sths.obj" \
	".\Release\auxi_cv.obj" \
	".\Release\comp_can.obj" \
	".\Release\restr_dp.obj" \
	".\Release\fix_recv.obj" \
	".\Release\sc_conls.obj" \
	".\Release\flux.obj" \
	".\Release\solu_com.obj" \
	".\Release\comp_dpn.obj" \
	".\Release\tools.obj" \
	".\Release\data_ext.obj" \
	".\Release\nd_ring.obj" \
	".\Release\isl_fix.obj" \
	".\Release\updcvcas_1n.obj" \
	".\Release\inputs_euler.obj" \
	".\Release\auxi_cc.obj" \
	".\Release\solution.obj" \
	".\Release\proc_pln.obj" \
	".\Release\comp_nd.obj" \
	".\Release\tools_rtn.obj" \
	".\Release\rmv-tgls.obj" \
	".\Release\pro_lcl.obj" \
	".\Release\auxi_cvp.obj" \
	".\Release\adss_cell.obj" \
	".\Release\reset.obj" \
	".\Release\grid_map.obj" \
	".\Release\solution_smooth.obj" \
	".\Release\flux_functions.obj" \
	".\Release\euler_functions.obj" \
	".\Release\smooth_xy.obj" \
	".\Release\pif_cell.obj" \
	".\Release\interpl.obj" \
	".\Release\fun_use.obj" \
	".\Release\output_for_show_euler.obj" \
	".\Release\upd_cvs.obj" \
	".\Release\pick_dp.obj" \
	".\Release\euler.obj" \
	".\Release\comp_on_cvn.obj" \
	".\Release\auxi_ccp.obj" \
	".\Release\updat_lr.obj" \
	".\Release\outputs_euler.obj" \
	".\Release\discontinuity_curve.obj" \
	".\Release\pro_smnd.obj" \
	".\Release\comp_in_smooth.obj" \
	".\Release\half_stp.obj" \
	".\Release\mesh_rt.obj" \
	".\Release\eva_sts.obj" \
	".\Release\updcvcas.obj" \
	".\Release\pick_mdp.obj" \
	".\Release\num_int.obj" \
	".\Release\data_pre.obj" \
	".\Release\dis_jntn.obj" \
	".\Release\crit_cell.obj" \
	".\Release\crit_info.obj" \
	".\Release\form_dsv.obj" \
	".\Release\out_mdl_euler.obj" \
	".\Release\ploynml.obj" \
	".\Release\num_int2.obj" \
	".\Release\nd_cls.obj" \
	".\Release\rcv_mpn.obj" \
	".\Release\locl_fix.obj" \
	".\Release\grid.obj" \
	".\Release\gif_cell.obj" \
	".\Release\bd_posi.obj" \
	".\Release\fif_cell.obj" \
	".\Release\sh_sts_el_HR.obj"

"$(OUTDIR)\main_euler.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\main_euler.exe" "$(OUTDIR)\critical_cell_information.mod"

CLEAN : 
	-@erase ".\Debug\critical_cell_information.mod"
	-@erase ".\Debug\adss_info_cell.mod"
	-@erase ".\Debug\geo_info_cell.mod"
	-@erase ".\Debug\grid.mod"
	-@erase ".\Debug\tools.mod"
	-@erase ".\Debug\polynomials.mod"
	-@erase ".\Debug\phy_info_cell.mod"
	-@erase ".\Debug\physical_state.mod"
	-@erase ".\Debug\euler_functions.mod"
	-@erase ".\Debug\data_extrapolation.mod"
	-@erase ".\Debug\main_euler.exe"
	-@erase ".\Debug\grid.obj"
	-@erase ".\Debug\euler_functions.obj"
	-@erase ".\Debug\comp_nd.obj"
	-@erase ".\Debug\node_cells.mod"
	-@erase ".\Debug\node_plug_ring.mod"
	-@erase ".\Debug\discontinuity_curves.mod"
	-@erase ".\Debug\critical_cells.mod"
	-@erase ".\Debug\auxiliary_discontinuity_curves.mod"
	-@erase ".\Debug\auxiliary_critical_cells.mod"
	-@erase ".\Debug\computational_information.mod"
	-@erase ".\Debug\grid_map.mod"
	-@erase ".\Debug\solu_in_smth.mod"
	-@erase ".\Debug\show_smth.mod"
	-@erase ".\Debug\pick_mdp.obj"
	-@erase ".\Debug\mesh_rt.obj"
	-@erase ".\Debug/solution.mod"
	-@erase ".\Debug\data_pre.obj"
	-@erase ".\Debug\solu_comput.mod"
	-@erase ".\Debug\dis_jntn.obj"
	-@erase ".\Debug\form_dsv.obj"
	-@erase ".\Debug\num_int2.obj"
	-@erase ".\Debug\numerical_integrals.mod"
	-@erase ".\Debug\function_in_use.mod"
	-@erase ".\Debug\crit_cell.obj"
	-@erase ".\Debug\gif_cell.obj"
	-@erase ".\Debug\ploynml.obj"
	-@erase ".\Debug\crit_info.obj"
	-@erase ".\Debug\flux.obj"
	-@erase ".\Debug\flux_functions.mod"
	-@erase ".\Debug\fif_cell.obj"
	-@erase ".\Debug\bd_posi.obj"
	-@erase ".\Debug\comp_can.obj"
	-@erase ".\Debug\restr_dp.obj"
	-@erase ".\Debug\main_temp.obj"
	-@erase ".\Debug/input_setup.mod"
	-@erase ".\Debug/out_intermediate.mod"
	-@erase ".\Debug/output_show.mod"
	-@erase ".\Debug\scan_conservation_on_list.mod"
	-@erase ".\Debug/auxiliary_curves_production.mod"
	-@erase ".\Debug\auxi_cc_prod.mod"
	-@erase ".\Debug\smoothen_xy_cell_averages.mod"
	-@erase ".\Debug\auxiliary_list_fix.mod"
	-@erase ".\Debug/computation_in_smooth.mod"
	-@erase ".\Debug\numerical_fluxes.mod"
	-@erase ".\Debug/computation_on_curves.mod"
	-@erase ".\Debug\compute_cell_average.mod"
	-@erase ".\Debug\recover_mid_positions.mod"
	-@erase ".\Debug\pick_discontinuity_positions.mod"
	-@erase ".\Debug\compute_discontinuity_positions.mod"
	-@erase ".\Debug\pick_middle_positions.mod"
	-@erase ".\Debug\restriction_of_positions.mod"
	-@erase ".\Debug\half_step_computation.mod"
	-@erase ".\Debug\compute_node_ca.mod"
	-@erase ".\Debug\output_middle_positions.mod"
	-@erase ".\Debug\compute_boundary_positions.mod"
	-@erase ".\Debug/reset_discvs.mod"
	-@erase ".\Debug\update_left_n_right_states.mod"
	-@erase ".\Debug\evaluate_updated_xxyy.mod"
	-@erase ".\Debug\produce_polynomial_4_xx_yy.mod"
	-@erase ".\Debug\interpolation.mod"
	-@erase ".\Debug\produce_local_smooth_4_nodes.mod"
	-@erase ".\Debug\produce_local_smooth.mod"
	-@erase ".\Debug\local_smooth_fix.mod"
	-@erase ".\Debug\data_preparation.mod"
	-@erase ".\Debug\update_curves.mod"
	-@erase ".\Debug\upd_cv_cases_1.mod"
	-@erase ".\Debug\upd_cv_cases.mod"
	-@erase ".\Debug\numerical_integral_in_2d.mod"
	-@erase ".\Debug\tools_for_reset.mod"
	-@erase ".\Debug\remove_tangles.mod"
	-@erase ".\Debug\fix_reversed_curves.mod"
	-@erase ".\Debug\form_discurves.mod"
	-@erase ".\Debug\disjoints_fix.mod"
	-@erase ".\Debug/mesh_ratio.mod"
	-@erase ".\Debug/output_show_HR.mod"
	-@erase ".\Debug\fix_recv.obj"
	-@erase ".\Debug\sh_sths.obj"
	-@erase ".\Debug\sc_conls.obj"
	-@erase ".\Debug\fun_use.obj"
	-@erase ".\Debug\euler.obj"
	-@erase ".\Debug\auxi_cv.obj"
	-@erase ".\Debug\data_ext.obj"
	-@erase ".\Debug\tools.obj"
	-@erase ".\Debug\updcvcas.obj"
	-@erase ".\Debug\solution.obj"
	-@erase ".\Debug\proc_pln.obj"
	-@erase ".\Debug\eva_sts.obj"
	-@erase ".\Debug\outputs_euler.obj"
	-@erase ".\Debug\auxi_cc.obj"
	-@erase ".\Debug\rmv-tgls.obj"
	-@erase ".\Debug\auxi_cvp.obj"
	-@erase ".\Debug\num_int.obj"
	-@erase ".\Debug\locl_fix.obj"
	-@erase ".\Debug\output_for_show_euler.obj"
	-@erase ".\Debug\pro_lcl.obj"
	-@erase ".\Debug\tools_rtn.obj"
	-@erase ".\Debug\rcv_mpn.obj"
	-@erase ".\Debug\adss_cell.obj"
	-@erase ".\Debug\reset.obj"
	-@erase ".\Debug\grid_map.obj"
	-@erase ".\Debug\out_mdl_euler.obj"
	-@erase ".\Debug\comp_in_smooth.obj"
	-@erase ".\Debug\nd_cls.obj"
	-@erase ".\Debug\comp_on_cvn.obj"
	-@erase ".\Debug\updcvcas_1n.obj"
	-@erase ".\Debug\solu_com.obj"
	-@erase ".\Debug\pif_cell.obj"
	-@erase ".\Debug\comp_dpn.obj"
	-@erase ".\Debug\discontinuity_curve.obj"
	-@erase ".\Debug\smooth_xy.obj"
	-@erase ".\Debug\interpl.obj"
	-@erase ".\Debug\solution_smooth.obj"
	-@erase ".\Debug\auxi_ccp.obj"
	-@erase ".\Debug\updat_lr.obj"
	-@erase ".\Debug\nd_ring.obj"
	-@erase ".\Debug\pro_smnd.obj"
	-@erase ".\Debug\inputs_euler.obj"
	-@erase ".\Debug\upd_cvs.obj"
	-@erase ".\Debug\pick_dp.obj"
	-@erase ".\Debug\isl_fix.obj"
	-@erase ".\Debug\half_stp.obj"
	-@erase ".\Debug\flux_functions.obj"
	-@erase ".\Debug\sh_sts_el_HR.obj"
	-@erase ".\Debug\main_euler.ilk"
	-@erase ".\Debug\main_euler.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Zi /I "Debug/" /c /nologo
F90_PROJ=/Zi /I "Debug/" /c /nologo /Fo"Debug/" /Fd"Debug/main_euler.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/main_euler.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/main_euler.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/main_euler.exe" 
LINK32_OBJS= \
	".\Debug\grid.obj" \
	".\Debug\euler_functions.obj" \
	".\Debug\comp_nd.obj" \
	".\Debug\pick_mdp.obj" \
	".\Debug\mesh_rt.obj" \
	".\Debug\data_pre.obj" \
	".\Debug\dis_jntn.obj" \
	".\Debug\form_dsv.obj" \
	".\Debug\num_int2.obj" \
	".\Debug\crit_cell.obj" \
	".\Debug\gif_cell.obj" \
	".\Debug\ploynml.obj" \
	".\Debug\crit_info.obj" \
	".\Debug\flux.obj" \
	".\Debug\fif_cell.obj" \
	".\Debug\bd_posi.obj" \
	".\Debug\comp_can.obj" \
	".\Debug\restr_dp.obj" \
	".\Debug\main_temp.obj" \
	".\Debug\fix_recv.obj" \
	".\Debug\sh_sths.obj" \
	".\Debug\sc_conls.obj" \
	".\Debug\fun_use.obj" \
	".\Debug\euler.obj" \
	".\Debug\auxi_cv.obj" \
	".\Debug\data_ext.obj" \
	".\Debug\tools.obj" \
	".\Debug\updcvcas.obj" \
	".\Debug\solution.obj" \
	".\Debug\proc_pln.obj" \
	".\Debug\eva_sts.obj" \
	".\Debug\outputs_euler.obj" \
	".\Debug\auxi_cc.obj" \
	".\Debug\rmv-tgls.obj" \
	".\Debug\auxi_cvp.obj" \
	".\Debug\num_int.obj" \
	".\Debug\locl_fix.obj" \
	".\Debug\output_for_show_euler.obj" \
	".\Debug\pro_lcl.obj" \
	".\Debug\tools_rtn.obj" \
	".\Debug\rcv_mpn.obj" \
	".\Debug\adss_cell.obj" \
	".\Debug\reset.obj" \
	".\Debug\grid_map.obj" \
	".\Debug\out_mdl_euler.obj" \
	".\Debug\comp_in_smooth.obj" \
	".\Debug\nd_cls.obj" \
	".\Debug\comp_on_cvn.obj" \
	".\Debug\updcvcas_1n.obj" \
	".\Debug\solu_com.obj" \
	".\Debug\pif_cell.obj" \
	".\Debug\comp_dpn.obj" \
	".\Debug\discontinuity_curve.obj" \
	".\Debug\smooth_xy.obj" \
	".\Debug\interpl.obj" \
	".\Debug\solution_smooth.obj" \
	".\Debug\auxi_ccp.obj" \
	".\Debug\updat_lr.obj" \
	".\Debug\nd_ring.obj" \
	".\Debug\pro_smnd.obj" \
	".\Debug\inputs_euler.obj" \
	".\Debug\upd_cvs.obj" \
	".\Debug\pick_dp.obj" \
	".\Debug\isl_fix.obj" \
	".\Debug\half_stp.obj" \
	".\Debug\flux_functions.obj" \
	".\Debug\sh_sts_el_HR.obj"

"$(OUTDIR)\main_euler.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "main_euler - Win32 Release"
# Name "main_euler - Win32 Debug"

!IF  "$(CFG)" == "main_euler - Win32 Release"

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\tools.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_TOOLS=\
	".\Release\polynomials.mod"\
	
F90_MODOUT=\
	"tools"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\tools.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\polynomials.mod"
   $(BuildCmds)

"$(INTDIR)\tools.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\polynomials.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_TOOLS=\
	".\Debug\polynomials.mod"\
	
F90_MODOUT=\
	"tools"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\tools.obj" : $(SOURCE) $(DEP_F90_TOOLS) "$(INTDIR)"\
 "$(INTDIR)\polynomials.mod"
   $(BuildCmds)

"$(INTDIR)\tools.mod" : $(SOURCE) $(DEP_F90_TOOLS) "$(INTDIR)"\
 "$(INTDIR)\polynomials.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\grid.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

F90_MODOUT=\
	"grid"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\grid.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

F90_MODOUT=\
	"grid"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\grid.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\euler\euler_functions.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

F90_MODOUT=\
	"euler_functions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\euler_functions.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\euler_functions.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

F90_MODOUT=\
	"euler_functions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\euler_functions.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\euler_functions.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\euler\euler.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_EULER=\
	".\Release\euler_functions.mod"\
	".\Release\data_extrapolation.mod"\
	
F90_MODOUT=\
	"physical_state"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\euler.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\euler_functions.mod"\
 "$(INTDIR)\data_extrapolation.mod"
   $(BuildCmds)

"$(INTDIR)\physical_state.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\euler_functions.mod" "$(INTDIR)\data_extrapolation.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_EULER=\
	".\Debug\euler_functions.mod"\
	".\Debug\data_extrapolation.mod"\
	
F90_MODOUT=\
	"physical_state"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\euler.obj" : $(SOURCE) $(DEP_F90_EULER) "$(INTDIR)"\
 "$(INTDIR)\euler_functions.mod" "$(INTDIR)\data_extrapolation.mod"
   $(BuildCmds)

"$(INTDIR)\physical_state.mod" : $(SOURCE) $(DEP_F90_EULER) "$(INTDIR)"\
 "$(INTDIR)\euler_functions.mod" "$(INTDIR)\data_extrapolation.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\solu\data_ext.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_DATA_=\
	".\Release\grid.mod"\
	
F90_MODOUT=\
	"data_extrapolation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\data_ext.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"
   $(BuildCmds)

"$(INTDIR)\data_extrapolation.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_DATA_=\
	".\Debug\grid.mod"\
	
F90_MODOUT=\
	"data_extrapolation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\data_ext.obj" : $(SOURCE) $(DEP_F90_DATA_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod"
   $(BuildCmds)

"$(INTDIR)\data_extrapolation.mod" : $(SOURCE) $(DEP_F90_DATA_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\fun_use.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_FUN_U=\
	".\Release\physical_state.mod"\
	
F90_MODOUT=\
	"function_in_use"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\fun_use.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

"$(INTDIR)\function_in_use.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_FUN_U=\
	".\Debug\physical_state.mod"\
	
F90_MODOUT=\
	"function_in_use"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\fun_use.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

"$(INTDIR)\function_in_use.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\adss_cell.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

F90_MODOUT=\
	"adss_info_cell"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\adss_cell.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\adss_info_cell.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

F90_MODOUT=\
	"adss_info_cell"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\adss_cell.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\adss_info_cell.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\gif_cell.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_GIF_C=\
	".\Release\grid.mod"\
	".\Release\tools.mod"\
	
F90_MODOUT=\
	"geo_info_cell"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\gif_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\tools.mod"
   $(BuildCmds)

"$(INTDIR)\geo_info_cell.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\tools.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_GIF_C=\
	".\Debug\grid.mod"\
	".\Debug\tools.mod"\
	
F90_MODOUT=\
	"geo_info_cell"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\gif_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\tools.mod"
   $(BuildCmds)

"$(INTDIR)\geo_info_cell.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\tools.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\pif_cell.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_PIF_C=\
	".\Release\physical_state.mod"\
	
F90_MODOUT=\
	"phy_info_cell"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pif_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

"$(INTDIR)\phy_info_cell.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_PIF_C=\
	".\Debug\physical_state.mod"\
	
F90_MODOUT=\
	"phy_info_cell"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pif_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

"$(INTDIR)\phy_info_cell.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\crit_info.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_CRIT_=\
	".\Release\adss_info_cell.mod"\
	".\Release\geo_info_cell.mod"\
	".\Release\phy_info_cell.mod"\
	
F90_MODOUT=\
	"critical_cell_information"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\crit_info.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\critical_cell_information.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_CRIT_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
F90_MODOUT=\
	"critical_cell_information"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\crit_info.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\critical_cell_information.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\reg_cuv\crit_cell.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_CRIT_C=\
	".\Release\geo_info_cell.mod"\
	".\Release\phy_info_cell.mod"\
	".\Release\adss_info_cell.mod"\
	
F90_MODOUT=\
	"critical_cells"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\crit_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod" "$(INTDIR)\adss_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\critical_cells.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\phy_info_cell.mod"\
 "$(INTDIR)\adss_info_cell.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_CRIT_C=\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	".\Debug\adss_info_cell.mod"\
	
F90_MODOUT=\
	"critical_cells"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\crit_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod" "$(INTDIR)\adss_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\critical_cells.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\phy_info_cell.mod"\
 "$(INTDIR)\adss_info_cell.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\reg_cuv\discontinuity_curve.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_DISCO=\
	".\Release\grid.mod"\
	".\Release\critical_cells.mod"\
	
F90_MODOUT=\
	"discontinuity_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\discontinuity_curve.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

"$(INTDIR)\discontinuity_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_DISCO=\
	".\Debug\grid.mod"\
	".\Debug\critical_cells.mod"\
	
F90_MODOUT=\
	"discontinuity_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\discontinuity_curve.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

"$(INTDIR)\discontinuity_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\auxi_cuv\fif_cell.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_FIF_C=\
	".\Release\physical_state.mod"\
	".\Release\grid.mod"\
	
F90_MODOUT=\
	"computational_information"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\fif_cell.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"\
 "$(INTDIR)\grid.mod"
   $(BuildCmds)

"$(INTDIR)\computational_information.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_FIF_C=\
	".\Debug\physical_state.mod"\
	".\Debug\grid.mod"\
	
F90_MODOUT=\
	"computational_information"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\fif_cell.obj" : $(SOURCE) $(DEP_F90_FIF_C) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod"
   $(BuildCmds)

"$(INTDIR)\computational_information.mod" : $(SOURCE) $(DEP_F90_FIF_C)\
 "$(INTDIR)" "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\auxi_cuv\auxi_cc.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_AUXI_=\
	".\Release\adss_info_cell.mod"\
	".\Release\geo_info_cell.mod"\
	".\Release\phy_info_cell.mod"\
	".\Release\computational_information.mod"\
	".\Release\critical_cells.mod"\
	
F90_MODOUT=\
	"auxiliary_critical_cells"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_cc.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\adss_info_cell.mod"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\phy_info_cell.mod"\
 "$(INTDIR)\computational_information.mod" "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_critical_cells.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod" "$(INTDIR)\computational_information.mod"\
 "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_AUXI_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	".\Debug\computational_information.mod"\
	".\Debug\critical_cells.mod"\
	
F90_MODOUT=\
	"auxiliary_critical_cells"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_cc.obj" : $(SOURCE) $(DEP_F90_AUXI_) "$(INTDIR)"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod" "$(INTDIR)\computational_information.mod"\
 "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_critical_cells.mod" : $(SOURCE) $(DEP_F90_AUXI_)\
 "$(INTDIR)" "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod" "$(INTDIR)\computational_information.mod"\
 "$(INTDIR)\critical_cells.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\dis_cuv\auxi_cuv\auxi_cv.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_AUXI_C=\
	".\Release\grid.mod"\
	".\Release\auxiliary_critical_cells.mod"\
	
F90_MODOUT=\
	"auxiliary_discontinuity_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_cv.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\auxiliary_critical_cells.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_discontinuity_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\auxiliary_critical_cells.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_AUXI_C=\
	".\Debug\grid.mod"\
	".\Debug\auxiliary_critical_cells.mod"\
	
F90_MODOUT=\
	"auxiliary_discontinuity_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_cv.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\auxiliary_critical_cells.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_discontinuity_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\auxiliary_critical_cells.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\euler\sh_sths.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SH_ST=\
	".\Release\physical_state.mod"\
	".\Release\euler_functions.mod"\
	
F90_MODOUT=\
	"show_smth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\sh_sths.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"\
 "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

"$(INTDIR)\show_smth.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SH_ST=\
	".\Debug\physical_state.mod"\
	".\Debug\euler_functions.mod"\
	
F90_MODOUT=\
	"show_smth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\sh_sths.obj" : $(SOURCE) $(DEP_F90_SH_ST) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

"$(INTDIR)\show_smth.mod" : $(SOURCE) $(DEP_F90_SH_ST) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\sl_smth\solution_smooth.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SOLUT=\
	".\Release\physical_state.mod"\
	".\Release\grid.mod"\
	".\Release\show_smth.mod"\
	
F90_MODOUT=\
	"solu_in_smth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solution_smooth.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod" "$(INTDIR)\show_smth.mod"
   $(BuildCmds)

"$(INTDIR)\solu_in_smth.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod" "$(INTDIR)\show_smth.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SOLUT=\
	".\Debug\physical_state.mod"\
	".\Debug\grid.mod"\
	".\Debug\show_smth.mod"\
	
F90_MODOUT=\
	"solu_in_smth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solution_smooth.obj" : $(SOURCE) $(DEP_F90_SOLUT) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod" "$(INTDIR)\show_smth.mod"
   $(BuildCmds)

"$(INTDIR)\solu_in_smth.mod" : $(SOURCE) $(DEP_F90_SOLUT) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod" "$(INTDIR)\show_smth.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\grid_map.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_GRID_=\
	".\Release\grid.mod"\
	".\Release\adss_info_cell.mod"\
	".\Release\discontinuity_curves.mod"\
	".\Release\auxiliary_discontinuity_curves.mod"\
	".\Release\solu_in_smth.mod"\
	
F90_MODOUT=\
	"grid_map"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\grid_map.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\solu_in_smth.mod"
   $(BuildCmds)

"$(INTDIR)\grid_map.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\adss_info_cell.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\solu_in_smth.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_GRID_=\
	".\Debug\grid.mod"\
	".\Debug\adss_info_cell.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\solu_in_smth.mod"\
	
F90_MODOUT=\
	"grid_map"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\grid_map.obj" : $(SOURCE) $(DEP_F90_GRID_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\adss_info_cell.mod"\
 "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\solu_in_smth.mod"
   $(BuildCmds)

"$(INTDIR)\grid_map.mod" : $(SOURCE) $(DEP_F90_GRID_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\adss_info_cell.mod"\
 "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\solu_in_smth.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\node\nd_ring.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_ND_RI=\
	".\Release\physical_state.mod"\
	".\Release\grid.mod"\
	
F90_MODOUT=\
	"node_plug_ring"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\nd_ring.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"\
 "$(INTDIR)\grid.mod"
   $(BuildCmds)

"$(INTDIR)\node_plug_ring.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_ND_RI=\
	".\Debug\physical_state.mod"\
	".\Debug\grid.mod"\
	
F90_MODOUT=\
	"node_plug_ring"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\nd_ring.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"\
 "$(INTDIR)\grid.mod"
   $(BuildCmds)

"$(INTDIR)\node_plug_ring.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\grid.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\node\nd_cls.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_ND_CL=\
	".\Release\node_plug_ring.mod"\
	".\Release\discontinuity_curves.mod"\
	
F90_MODOUT=\
	"node_cells"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\nd_cls.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\node_plug_ring.mod"\
 "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\node_cells.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\node_plug_ring.mod" "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_ND_CL=\
	".\Debug\node_plug_ring.mod"\
	".\Debug\discontinuity_curves.mod"\
	
F90_MODOUT=\
	"node_cells"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\nd_cls.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\node_plug_ring.mod"\
 "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\node_cells.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\node_plug_ring.mod" "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\solu\solution.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SOLUTI=\
	".\Release\grid_map.mod"\
	".\Release\solu_in_smth.mod"\
	".\Release\discontinuity_curves.mod"\
	".\Release\node_cells.mod"\
	
F90_MODOUT=\
	"solution"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solution.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid_map.mod"\
 "$(INTDIR)\solu_in_smth.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

"$(INTDIR)\solution.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid_map.mod"\
 "$(INTDIR)\solu_in_smth.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SOLUTI=\
	".\Debug\grid_map.mod"\
	".\Debug\solu_in_smth.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\node_cells.mod"\
	
F90_MODOUT=\
	"solution"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solution.obj" : $(SOURCE) $(DEP_F90_SOLUTI) "$(INTDIR)"\
 "$(INTDIR)\grid_map.mod" "$(INTDIR)\solu_in_smth.mod"\
 "$(INTDIR)\discontinuity_curves.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

"$(INTDIR)\solution.mod" : $(SOURCE) $(DEP_F90_SOLUTI) "$(INTDIR)"\
 "$(INTDIR)\grid_map.mod" "$(INTDIR)\solu_in_smth.mod"\
 "$(INTDIR)\discontinuity_curves.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\interface\outputs_euler.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_OUTPU=\
	".\Release\solution.mod"\
	".\Release\grid_map.mod"\
	
F90_MODOUT=\
	"out_intermediate"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\outputs_euler.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"\
 "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

"$(INTDIR)\out_intermediate.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_OUTPU=\
	".\Debug/solution.mod"\
	".\Debug\grid_map.mod"\
	
F90_MODOUT=\
	"out_intermediate"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\outputs_euler.obj" : $(SOURCE) $(DEP_F90_OUTPU) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

"$(INTDIR)\out_intermediate.mod" : $(SOURCE) $(DEP_F90_OUTPU) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\main_itf\main_temp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_MAIN_=\
	".\Release\input_setup.mod"\
	".\Release\out_intermediate.mod"\
	".\Release\output_show.mod"\
	".\Release\auxiliary_curves_production.mod"\
	".\Release\computation_in_smooth.mod"\
	".\Release\computation_on_curves.mod"\
	".\Release\reset_discvs.mod"\
	".\Release\mesh_ratio.mod"\
	".\Release\output_show_HR.mod"\
	

"$(INTDIR)\main_temp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\input_setup.mod"\
 "$(INTDIR)\out_intermediate.mod" "$(INTDIR)\output_show.mod"\
 "$(INTDIR)\auxiliary_curves_production.mod"\
 "$(INTDIR)\computation_in_smooth.mod" "$(INTDIR)\computation_on_curves.mod"\
 "$(INTDIR)\reset_discvs.mod" "$(INTDIR)\mesh_ratio.mod"\
 "$(INTDIR)\output_show_HR.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_MAIN_=\
	".\Debug/input_setup.mod"\
	".\Debug/out_intermediate.mod"\
	".\Debug/output_show.mod"\
	".\Debug/auxiliary_curves_production.mod"\
	".\Debug/computation_in_smooth.mod"\
	".\Debug/computation_on_curves.mod"\
	".\Debug/reset_discvs.mod"\
	".\Debug/mesh_ratio.mod"\
	".\Debug/output_show_HR.mod"\
	

"$(INTDIR)\main_temp.obj" : $(SOURCE) $(DEP_F90_MAIN_) "$(INTDIR)"\
 "$(INTDIR)\input_setup.mod" "$(INTDIR)\out_intermediate.mod"\
 "$(INTDIR)\output_show.mod" "$(INTDIR)\auxiliary_curves_production.mod"\
 "$(INTDIR)\computation_in_smooth.mod" "$(INTDIR)\computation_on_curves.mod"\
 "$(INTDIR)\reset_discvs.mod" "$(INTDIR)\mesh_ratio.mod"\
 "$(INTDIR)\output_show_HR.mod"
   $(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\interface\inputs_euler.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_INPUT=\
	".\Release\solution.mod"\
	
F90_MODOUT=\
	"input_setup"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\inputs_euler.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\input_setup.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_INPUT=\
	".\Debug/solution.mod"\
	
F90_MODOUT=\
	"input_setup"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\inputs_euler.obj" : $(SOURCE) $(DEP_F90_INPUT) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\input_setup.mod" : $(SOURCE) $(DEP_F90_INPUT) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\interface\output_for_show_euler.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_OUTPUT=\
	".\Release\solution.mod"\
	".\Release\euler_functions.mod"\
	".\Release\scan_conservation_on_list.mod"\
	
F90_MODOUT=\
	"output_show"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\output_for_show_euler.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\euler_functions.mod"\
 "$(INTDIR)\scan_conservation_on_list.mod"
   $(BuildCmds)

"$(INTDIR)\output_show.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"\
 "$(INTDIR)\euler_functions.mod" "$(INTDIR)\scan_conservation_on_list.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_OUTPUT=\
	".\Debug/solution.mod"\
	".\Debug\euler_functions.mod"\
	".\Debug\scan_conservation_on_list.mod"\
	
F90_MODOUT=\
	"output_show"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\output_for_show_euler.obj" : $(SOURCE) $(DEP_F90_OUTPUT) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\euler_functions.mod"\
 "$(INTDIR)\scan_conservation_on_list.mod"
   $(BuildCmds)

"$(INTDIR)\output_show.mod" : $(SOURCE) $(DEP_F90_OUTPUT) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\euler_functions.mod"\
 "$(INTDIR)\scan_conservation_on_list.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\mesh_ratio\mesh_rt.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_MESH_=\
	".\Release\solution.mod"\
	
F90_MODOUT=\
	"mesh_ratio"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\mesh_rt.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\mesh_ratio.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_MESH_=\
	".\Debug/solution.mod"\
	
F90_MODOUT=\
	"mesh_ratio"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\mesh_rt.obj" : $(SOURCE) $(DEP_F90_MESH_) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\mesh_ratio.mod" : $(SOURCE) $(DEP_F90_MESH_) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\prc_auxi\auxi_ccp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_AUXI_CC=\
	".\Release\geo_info_cell.mod"\
	".\Release\phy_info_cell.mod"\
	
F90_MODOUT=\
	"auxi_cc_prod"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_ccp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\auxi_cc_prod.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_AUXI_CC=\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
F90_MODOUT=\
	"auxi_cc_prod"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_ccp.obj" : $(SOURCE) $(DEP_F90_AUXI_CC) "$(INTDIR)"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\auxi_cc_prod.mod" : $(SOURCE) $(DEP_F90_AUXI_CC) "$(INTDIR)"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\phy_info_cell.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\smooth_xy.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SMOOT=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"smoothen_xy_cell_averages"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\smooth_xy.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\smoothen_xy_cell_averages.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SMOOT=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"smoothen_xy_cell_averages"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\smooth_xy.obj" : $(SOURCE) $(DEP_F90_SMOOT) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\smoothen_xy_cell_averages.mod" : $(SOURCE) $(DEP_F90_SMOOT)\
 "$(INTDIR)" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\prc_auxi\auxi_cvp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_AUXI_CV=\
	".\Release\discontinuity_curves.mod"\
	".\Release\auxiliary_discontinuity_curves.mod"\
	".\Release\auxi_cc_prod.mod"\
	".\Release\auxiliary_list_fix.mod"\
	".\Release\node_cells.mod"\
	".\Release\smoothen_xy_cell_averages.mod"\
	
F90_MODOUT=\
	"auxiliary_curves_production"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_cvp.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\smoothen_xy_cell_averages.mod"\
 "$(INTDIR)\auxiliary_list_fix.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_curves_production.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\smoothen_xy_cell_averages.mod"\
 "$(INTDIR)\auxiliary_list_fix.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_AUXI_CV=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\auxi_cc_prod.mod"\
	".\Debug\auxiliary_list_fix.mod"\
	".\Debug\node_cells.mod"\
	".\Debug\smoothen_xy_cell_averages.mod"\
	
F90_MODOUT=\
	"auxiliary_curves_production"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\auxi_cvp.obj" : $(SOURCE) $(DEP_F90_AUXI_CV) "$(INTDIR)"\
 "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\smoothen_xy_cell_averages.mod"\
 "$(INTDIR)\auxiliary_list_fix.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_curves_production.mod" : $(SOURCE) $(DEP_F90_AUXI_CV)\
 "$(INTDIR)" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\smoothen_xy_cell_averages.mod"\
 "$(INTDIR)\auxiliary_list_fix.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\phys_flx\flux_functions.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_FLUX_=\
	".\Release\grid.mod"\
	".\Release\euler_functions.mod"\
	
F90_MODOUT=\
	"flux_functions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\flux_functions.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

"$(INTDIR)\flux_functions.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"\
 "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_FLUX_=\
	".\Debug\grid.mod"\
	".\Debug\euler_functions.mod"\
	
F90_MODOUT=\
	"flux_functions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\flux_functions.obj" : $(SOURCE) $(DEP_F90_FLUX_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

"$(INTDIR)\flux_functions.mod" : $(SOURCE) $(DEP_F90_FLUX_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\phys_flx\flux.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_FLUX_F=\
	".\Release\physical_state.mod"\
	".\Release\flux_functions.mod"\
	
F90_MODOUT=\
	"numerical_fluxes"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\flux.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"\
 "$(INTDIR)\flux_functions.mod"
   $(BuildCmds)

"$(INTDIR)\numerical_fluxes.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\flux_functions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_FLUX_F=\
	".\Debug\physical_state.mod"\
	".\Debug\flux_functions.mod"\
	
F90_MODOUT=\
	"numerical_fluxes"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\flux.obj" : $(SOURCE) $(DEP_F90_FLUX_F) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\flux_functions.mod"
   $(BuildCmds)

"$(INTDIR)\numerical_fluxes.mod" : $(SOURCE) $(DEP_F90_FLUX_F) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\flux_functions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_sth\comp_in_smooth.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_COMP_=\
	".\Release\solu_comput.mod"\
	".\Release\numerical_fluxes.mod"\
	".\Release\node_cells.mod"\
	
F90_MODOUT=\
	"computation_in_smooth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_in_smooth.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\numerical_fluxes.mod" "$(INTDIR)\node_cells.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\computation_in_smooth.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\numerical_fluxes.mod" "$(INTDIR)\node_cells.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_COMP_=\
	".\Debug\solu_comput.mod"\
	".\Debug\numerical_fluxes.mod"\
	".\Debug\node_cells.mod"\
	
F90_MODOUT=\
	"computation_in_smooth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_in_smooth.obj" : $(SOURCE) $(DEP_F90_COMP_) "$(INTDIR)"\
 "$(INTDIR)\numerical_fluxes.mod" "$(INTDIR)\node_cells.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\computation_in_smooth.mod" : $(SOURCE) $(DEP_F90_COMP_) "$(INTDIR)"\
 "$(INTDIR)\numerical_fluxes.mod" "$(INTDIR)\node_cells.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\solu\solu_com.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SOLU_=\
	".\Release\grid_map.mod"\
	".\Release\solu_in_smth.mod"\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"solu_comput"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solu_com.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid_map.mod"\
 "$(INTDIR)\solu_in_smth.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\solu_comput.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid_map.mod"\
 "$(INTDIR)\solu_in_smth.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SOLU_=\
	".\Debug\grid_map.mod"\
	".\Debug\solu_in_smth.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"solu_comput"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\solu_com.obj" : $(SOURCE) $(DEP_F90_SOLU_) "$(INTDIR)"\
 "$(INTDIR)\grid_map.mod" "$(INTDIR)\solu_in_smth.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\solu_comput.mod" : $(SOURCE) $(DEP_F90_SOLU_) "$(INTDIR)"\
 "$(INTDIR)\grid_map.mod" "$(INTDIR)\solu_in_smth.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\pick_dp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_PICK_=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"pick_discontinuity_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pick_dp.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\pick_discontinuity_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_PICK_=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"pick_discontinuity_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pick_dp.obj" : $(SOURCE) $(DEP_F90_PICK_) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\pick_discontinuity_positions.mod" : $(SOURCE) $(DEP_F90_PICK_)\
 "$(INTDIR)" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\rcv_mpn.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_RCV_M=\
	".\Release\solu_comput.mod"\
	".\Release\pick_discontinuity_positions.mod"\
	
F90_MODOUT=\
	"recover_mid_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\rcv_mpn.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\pick_discontinuity_positions.mod"
   $(BuildCmds)

"$(INTDIR)\recover_mid_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\pick_discontinuity_positions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_RCV_M=\
	".\Debug\solu_comput.mod"\
	".\Debug\pick_discontinuity_positions.mod"\
	
F90_MODOUT=\
	"recover_mid_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\rcv_mpn.obj" : $(SOURCE) $(DEP_F90_RCV_M) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\pick_discontinuity_positions.mod"
   $(BuildCmds)

"$(INTDIR)\recover_mid_positions.mod" : $(SOURCE) $(DEP_F90_RCV_M) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\pick_discontinuity_positions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\pick_mdp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_PICK_M=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"pick_middle_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pick_mdp.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\pick_middle_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_PICK_M=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"pick_middle_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pick_mdp.obj" : $(SOURCE) $(DEP_F90_PICK_M) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\pick_middle_positions.mod" : $(SOURCE) $(DEP_F90_PICK_M) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_dpn.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_COMP_D=\
	".\Release\pick_middle_positions.mod"\
	
F90_MODOUT=\
	"compute_discontinuity_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_dpn.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\pick_middle_positions.mod"
   $(BuildCmds)

"$(INTDIR)\compute_discontinuity_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\pick_middle_positions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_COMP_D=\
	".\Debug\pick_middle_positions.mod"\
	
F90_MODOUT=\
	"compute_discontinuity_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_dpn.obj" : $(SOURCE) $(DEP_F90_COMP_D) "$(INTDIR)"\
 "$(INTDIR)\pick_middle_positions.mod"
   $(BuildCmds)

"$(INTDIR)\compute_discontinuity_positions.mod" : $(SOURCE) $(DEP_F90_COMP_D)\
 "$(INTDIR)" "$(INTDIR)\pick_middle_positions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\restr_dp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_RESTR=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"restriction_of_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\restr_dp.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\restriction_of_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_RESTR=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"restriction_of_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\restr_dp.obj" : $(SOURCE) $(DEP_F90_RESTR) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\restriction_of_positions.mod" : $(SOURCE) $(DEP_F90_RESTR)\
 "$(INTDIR)" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\half_stp.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_HALF_=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"half_step_computation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\half_stp.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\half_step_computation.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_HALF_=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"half_step_computation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\half_stp.obj" : $(SOURCE) $(DEP_F90_HALF_) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\half_step_computation.mod" : $(SOURCE) $(DEP_F90_HALF_) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_can.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_COMP_C=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	".\Release\geo_info_cell.mod"\
	
F90_MODOUT=\
	"compute_cell_average"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_can.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\geo_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\compute_cell_average.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\geo_info_cell.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_COMP_C=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\geo_info_cell.mod"\
	
F90_MODOUT=\
	"compute_cell_average"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_can.obj" : $(SOURCE) $(DEP_F90_COMP_C) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\geo_info_cell.mod"
   $(BuildCmds)

"$(INTDIR)\compute_cell_average.mod" : $(SOURCE) $(DEP_F90_COMP_C) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\geo_info_cell.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_nd.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_COMP_N=\
	".\Release\node_cells.mod"\
	".\Release\auxiliary_discontinuity_curves.mod"\
	".\Release\geo_info_cell.mod"\
	".\Release\grid_map.mod"\
	
F90_MODOUT=\
	"compute_node_ca"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_nd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\node_cells.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\geo_info_cell.mod"\
 "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

"$(INTDIR)\compute_node_ca.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_COMP_N=\
	".\Debug\node_cells.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\grid_map.mod"\
	
F90_MODOUT=\
	"compute_node_ca"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_nd.obj" : $(SOURCE) $(DEP_F90_COMP_N) "$(INTDIR)"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

"$(INTDIR)\compute_node_ca.mod" : $(SOURCE) $(DEP_F90_COMP_N) "$(INTDIR)"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"\
 "$(INTDIR)\geo_info_cell.mod" "$(INTDIR)\grid_map.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\interface\out_mdl_euler.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_OUT_M=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"output_middle_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\out_mdl_euler.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\output_middle_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_OUT_M=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"output_middle_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\out_mdl_euler.obj" : $(SOURCE) $(DEP_F90_OUT_M) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\output_middle_positions.mod" : $(SOURCE) $(DEP_F90_OUT_M)\
 "$(INTDIR)" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\comp_on_cvn.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_COMP_O=\
	".\Release\compute_cell_average.mod"\
	".\Release\recover_mid_positions.mod"\
	".\Release\compute_discontinuity_positions.mod"\
	".\Release\compute_boundary_positions.mod"\
	".\Release\restriction_of_positions.mod"\
	".\Release\half_step_computation.mod"\
	".\Release\compute_node_ca.mod"\
	".\Release\output_middle_positions.mod"\
	
F90_MODOUT=\
	"computation_on_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_on_cvn.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\compute_cell_average.mod" "$(INTDIR)\recover_mid_positions.mod"\
 "$(INTDIR)\compute_discontinuity_positions.mod"\
 "$(INTDIR)\restriction_of_positions.mod" "$(INTDIR)\half_step_computation.mod"\
 "$(INTDIR)\compute_node_ca.mod" "$(INTDIR)\output_middle_positions.mod"\
 "$(INTDIR)\compute_boundary_positions.mod"
   $(BuildCmds)

"$(INTDIR)\computation_on_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\compute_cell_average.mod" "$(INTDIR)\recover_mid_positions.mod"\
 "$(INTDIR)\compute_discontinuity_positions.mod"\
 "$(INTDIR)\restriction_of_positions.mod" "$(INTDIR)\half_step_computation.mod"\
 "$(INTDIR)\compute_node_ca.mod" "$(INTDIR)\output_middle_positions.mod"\
 "$(INTDIR)\compute_boundary_positions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_COMP_O=\
	".\Debug\compute_cell_average.mod"\
	".\Debug\recover_mid_positions.mod"\
	".\Debug\compute_discontinuity_positions.mod"\
	".\Debug\compute_boundary_positions.mod"\
	".\Debug\restriction_of_positions.mod"\
	".\Debug\half_step_computation.mod"\
	".\Debug\compute_node_ca.mod"\
	".\Debug\output_middle_positions.mod"\
	
F90_MODOUT=\
	"computation_on_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\comp_on_cvn.obj" : $(SOURCE) $(DEP_F90_COMP_O) "$(INTDIR)"\
 "$(INTDIR)\compute_cell_average.mod" "$(INTDIR)\recover_mid_positions.mod"\
 "$(INTDIR)\compute_discontinuity_positions.mod"\
 "$(INTDIR)\restriction_of_positions.mod" "$(INTDIR)\half_step_computation.mod"\
 "$(INTDIR)\compute_node_ca.mod" "$(INTDIR)\output_middle_positions.mod"\
 "$(INTDIR)\compute_boundary_positions.mod"
   $(BuildCmds)

"$(INTDIR)\computation_on_curves.mod" : $(SOURCE) $(DEP_F90_COMP_O) "$(INTDIR)"\
 "$(INTDIR)\compute_cell_average.mod" "$(INTDIR)\recover_mid_positions.mod"\
 "$(INTDIR)\compute_discontinuity_positions.mod"\
 "$(INTDIR)\restriction_of_positions.mod" "$(INTDIR)\half_step_computation.mod"\
 "$(INTDIR)\compute_node_ca.mod" "$(INTDIR)\output_middle_positions.mod"\
 "$(INTDIR)\compute_boundary_positions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\comp_cv\bd_posi.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_BD_PO=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"compute_boundary_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\bd_posi.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\compute_boundary_positions.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_BD_PO=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"compute_boundary_positions"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\bd_posi.obj" : $(SOURCE) $(DEP_F90_BD_PO) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\compute_boundary_positions.mod" : $(SOURCE) $(DEP_F90_BD_PO)\
 "$(INTDIR)" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\tools\interpl.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

F90_MODOUT=\
	"interpolation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\interpl.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\interpolation.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

F90_MODOUT=\
	"interpolation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\interpl.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\interpolation.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\phys_flx\proc_pln.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_PROC_=\
	".\Release\interpolation.mod"\
	".\Release\solu_comput.mod"\
	
F90_MODOUT=\
	"produce_polynomial_4_xx_yy"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\proc_pln.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\interpolation.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\produce_polynomial_4_xx_yy.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\interpolation.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_PROC_=\
	".\Debug\interpolation.mod"\
	".\Debug\solu_comput.mod"\
	
F90_MODOUT=\
	"produce_polynomial_4_xx_yy"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\proc_pln.obj" : $(SOURCE) $(DEP_F90_PROC_) "$(INTDIR)"\
 "$(INTDIR)\interpolation.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\produce_polynomial_4_xx_yy.mod" : $(SOURCE) $(DEP_F90_PROC_)\
 "$(INTDIR)" "$(INTDIR)\interpolation.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\phys_flx\phy_updt\updat_lr.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_UPDAT=\
	".\Release\evaluate_updated_xxyy.mod"\
	
F90_MODOUT=\
	"update_left_n_right_states"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\updat_lr.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\evaluate_updated_xxyy.mod"
   $(BuildCmds)

"$(INTDIR)\update_left_n_right_states.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\evaluate_updated_xxyy.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_UPDAT=\
	".\Debug\evaluate_updated_xxyy.mod"\
	
F90_MODOUT=\
	"update_left_n_right_states"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\updat_lr.obj" : $(SOURCE) $(DEP_F90_UPDAT) "$(INTDIR)"\
 "$(INTDIR)\evaluate_updated_xxyy.mod"
   $(BuildCmds)

"$(INTDIR)\update_left_n_right_states.mod" : $(SOURCE) $(DEP_F90_UPDAT)\
 "$(INTDIR)" "$(INTDIR)\evaluate_updated_xxyy.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\interact\pro_lcl.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_PRO_L=\
	".\Release\solution.mod"\
	".\Release\solu_comput.mod"\
	
F90_MODOUT=\
	"produce_local_smooth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pro_lcl.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\produce_local_smooth.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_PRO_L=\
	".\Debug/solution.mod"\
	".\Debug\solu_comput.mod"\
	
F90_MODOUT=\
	"produce_local_smooth"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pro_lcl.obj" : $(SOURCE) $(DEP_F90_PRO_L) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\produce_local_smooth.mod" : $(SOURCE) $(DEP_F90_PRO_L) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\pro_smnd.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_PRO_S=\
	".\Release\produce_local_smooth.mod"\
	
F90_MODOUT=\
	"produce_local_smooth_4_nodes"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pro_smnd.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\produce_local_smooth.mod"
   $(BuildCmds)

"$(INTDIR)\produce_local_smooth_4_nodes.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\produce_local_smooth.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

NODEP_F90_PRO_S=\
	".\Debug\produce_local_smooth.mod"\
	
F90_MODOUT=\
	"produce_local_smooth_4_nodes"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\pro_smnd.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\produce_local_smooth.mod"
   $(BuildCmds)

"$(INTDIR)\produce_local_smooth_4_nodes.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\produce_local_smooth.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\locl_fix.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_LOCL_=\
	".\Release\solution.mod"\
	".\Release\solu_comput.mod"\
	
F90_MODOUT=\
	"local_smooth_fix"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\locl_fix.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"\
 "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\local_smooth_fix.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_LOCL_=\
	".\Debug/solution.mod"\
	".\Debug\solu_comput.mod"\
	
F90_MODOUT=\
	"local_smooth_fix"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\locl_fix.obj" : $(SOURCE) $(DEP_F90_LOCL_) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

"$(INTDIR)\local_smooth_fix.mod" : $(SOURCE) $(DEP_F90_LOCL_) "$(INTDIR)"\
 "$(INTDIR)\solution.mod" "$(INTDIR)\solu_comput.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\data_pre.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_DATA_P=\
	".\Release\solu_comput.mod"\
	".\Release\node_cells.mod"\
	
F90_MODOUT=\
	"data_preparation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\data_pre.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

"$(INTDIR)\data_preparation.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_DATA_P=\
	".\Debug\solu_comput.mod"\
	".\Debug\node_cells.mod"\
	
F90_MODOUT=\
	"data_preparation"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\data_pre.obj" : $(SOURCE) $(DEP_F90_DATA_P) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

"$(INTDIR)\data_preparation.mod" : $(SOURCE) $(DEP_F90_DATA_P) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\updcvcas_1n.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_UPDCV=\
	".\Release\upd_cv_cases.mod"\
	".\Release\node_cells.mod"\
	
F90_MODOUT=\
	"upd_cv_cases_1"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\updcvcas_1n.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\node_cells.mod"\
 "$(INTDIR)\upd_cv_cases.mod"
   $(BuildCmds)

"$(INTDIR)\upd_cv_cases_1.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\upd_cv_cases.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_UPDCV=\
	".\Debug\node_cells.mod"\
	
NODEP_F90_UPDCV=\
	".\Debug\upd_cv_cases.mod"\
	
F90_MODOUT=\
	"upd_cv_cases_1"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\updcvcas_1n.obj" : $(SOURCE) $(DEP_F90_UPDCV) "$(INTDIR)"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\upd_cv_cases.mod"
   $(BuildCmds)

"$(INTDIR)\upd_cv_cases_1.mod" : $(SOURCE) $(DEP_F90_UPDCV) "$(INTDIR)"\
 "$(INTDIR)\node_cells.mod" "$(INTDIR)\upd_cv_cases.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\updcvcas.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_UPDCVC=\
	".\Release\solu_comput.mod"\
	".\Release\discontinuity_curves.mod"\
	".\Release\tools_for_reset.mod"\
	".\Release\numerical_integral_in_2D.mod"\
	
F90_MODOUT=\
	"upd_cv_cases"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\updcvcas.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\discontinuity_curves.mod" "$(INTDIR)\numerical_integral_in_2D.mod"\
 "$(INTDIR)\tools_for_reset.mod"
   $(BuildCmds)

"$(INTDIR)\upd_cv_cases.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\numerical_integral_in_2D.mod" "$(INTDIR)\tools_for_reset.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_UPDCVC=\
	".\Debug\solu_comput.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\tools_for_reset.mod"\
	".\Debug\numerical_integral_in_2d.mod"\
	
F90_MODOUT=\
	"upd_cv_cases"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\updcvcas.obj" : $(SOURCE) $(DEP_F90_UPDCVC) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\numerical_integral_in_2d.mod" "$(INTDIR)\tools_for_reset.mod"
   $(BuildCmds)

"$(INTDIR)\upd_cv_cases.mod" : $(SOURCE) $(DEP_F90_UPDCVC) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\discontinuity_curves.mod"\
 "$(INTDIR)\numerical_integral_in_2d.mod" "$(INTDIR)\tools_for_reset.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE="\workf90_1\bs_oprt\rest_cv\rmv-tgls.f90"

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_RMV_T=\
	".\Release\solu_comput.mod"\
	".\Release\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"remove_tangles"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\rmv-tgls.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\remove_tangles.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_RMV_T=\
	".\Debug\solu_comput.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	
F90_MODOUT=\
	"remove_tangles"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\rmv-tgls.obj" : $(SOURCE) $(DEP_F90_RMV_T) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\remove_tangles.mod" : $(SOURCE) $(DEP_F90_RMV_T) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\auxiliary_discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\num_int.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_NUM_I=\
	".\Release\physical_state.mod"\
	
F90_MODOUT=\
	"numerical_integrals"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\num_int.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

"$(INTDIR)\numerical_integrals.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_NUM_I=\
	".\Debug\physical_state.mod"\
	
F90_MODOUT=\
	"numerical_integrals"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\num_int.obj" : $(SOURCE) $(DEP_F90_NUM_I) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

"$(INTDIR)\numerical_integrals.mod" : $(SOURCE) $(DEP_F90_NUM_I) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\num_int2.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_NUM_IN=\
	".\Release\physical_state.mod"\
	".\Release\numerical_integrals.mod"\
	".\Release\function_in_use.mod"\
	
F90_MODOUT=\
	"numerical_integral_in_2D"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\num_int2.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\physical_state.mod"\
 "$(INTDIR)\numerical_integrals.mod" "$(INTDIR)\function_in_use.mod"
   $(BuildCmds)

"$(INTDIR)\numerical_integral_in_2D.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\numerical_integrals.mod"\
 "$(INTDIR)\function_in_use.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_NUM_IN=\
	".\Debug\physical_state.mod"\
	".\Debug\numerical_integrals.mod"\
	".\Debug\function_in_use.mod"\
	
F90_MODOUT=\
	"numerical_integral_in_2D"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\num_int2.obj" : $(SOURCE) $(DEP_F90_NUM_IN) "$(INTDIR)"\
 "$(INTDIR)\physical_state.mod" "$(INTDIR)\numerical_integrals.mod"\
 "$(INTDIR)\function_in_use.mod"
   $(BuildCmds)

"$(INTDIR)\numerical_integral_in_2d.mod" : $(SOURCE) $(DEP_F90_NUM_IN)\
 "$(INTDIR)" "$(INTDIR)\physical_state.mod" "$(INTDIR)\numerical_integrals.mod"\
 "$(INTDIR)\function_in_use.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\ploynml.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

F90_MODOUT=\
	"polynomials"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\ploynml.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\polynomials.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

F90_MODOUT=\
	"polynomials"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\ploynml.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\polynomials.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\upd_cvs.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_UPD_C=\
	".\Release\upd_cv_cases_1.mod"\
	".\Release\remove_tangles.mod"\
	
F90_MODOUT=\
	"update_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\upd_cvs.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\upd_cv_cases_1.mod"\
 "$(INTDIR)\remove_tangles.mod"
   $(BuildCmds)

"$(INTDIR)\update_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\upd_cv_cases_1.mod" "$(INTDIR)\remove_tangles.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_UPD_C=\
	".\Debug\upd_cv_cases_1.mod"\
	".\Debug\remove_tangles.mod"\
	
F90_MODOUT=\
	"update_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\upd_cvs.obj" : $(SOURCE) $(DEP_F90_UPD_C) "$(INTDIR)"\
 "$(INTDIR)\upd_cv_cases_1.mod" "$(INTDIR)\remove_tangles.mod"
   $(BuildCmds)

"$(INTDIR)\update_curves.mod" : $(SOURCE) $(DEP_F90_UPD_C) "$(INTDIR)"\
 "$(INTDIR)\upd_cv_cases_1.mod" "$(INTDIR)\remove_tangles.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\fix_recv.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_FIX_R=\
	".\Release\solu_comput.mod"\
	".\Release\node_cells.mod"\
	
F90_MODOUT=\
	"fix_reversed_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\fix_recv.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

"$(INTDIR)\fix_reversed_curves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_FIX_R=\
	".\Debug\solu_comput.mod"\
	".\Debug\node_cells.mod"\
	
F90_MODOUT=\
	"fix_reversed_curves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\fix_recv.obj" : $(SOURCE) $(DEP_F90_FIX_R) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

"$(INTDIR)\fix_reversed_curves.mod" : $(SOURCE) $(DEP_F90_FIX_R) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\node_cells.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\tools_rtn.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_TOOLS_=\
	".\Release\solu_comput.mod"\
	".\Release\tools.mod"\
	
F90_MODOUT=\
	"tools_for_reset"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\tools_rtn.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\tools.mod"
   $(BuildCmds)

"$(INTDIR)\tools_for_reset.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\tools.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_TOOLS_=\
	".\Debug\solu_comput.mod"\
	".\Debug\tools.mod"\
	
F90_MODOUT=\
	"tools_for_reset"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\tools_rtn.obj" : $(SOURCE) $(DEP_F90_TOOLS_) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\tools.mod"
   $(BuildCmds)

"$(INTDIR)\tools_for_reset.mod" : $(SOURCE) $(DEP_F90_TOOLS_) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\tools.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\form_dsv.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_FORM_=\
	".\Release\solu_comput.mod"\
	".\Release\discontinuity_curves.mod"\
	
F90_MODOUT=\
	"form_discurves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\form_dsv.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solu_comput.mod"\
 "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\form_discurves.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_FORM_=\
	".\Debug\solu_comput.mod"\
	".\Debug\discontinuity_curves.mod"\
	
F90_MODOUT=\
	"form_discurves"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\form_dsv.obj" : $(SOURCE) $(DEP_F90_FORM_) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\form_discurves.mod" : $(SOURCE) $(DEP_F90_FORM_) "$(INTDIR)"\
 "$(INTDIR)\solu_comput.mod" "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\rest_cv\reset.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_RESET=\
	".\Release\update_left_n_right_states.mod"\
	".\Release\produce_local_smooth_4_nodes.mod"\
	".\Release\local_smooth_fix.mod"\
	".\Release\data_preparation.mod"\
	".\Release\update_curves.mod"\
	".\Release\fix_reversed_curves.mod"\
	".\Release\form_discurves.mod"\
	".\Release\disjoints_fix.mod"\
	
F90_MODOUT=\
	"reset_discvs"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\reset.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\update_left_n_right_states.mod"\
 "$(INTDIR)\produce_local_smooth_4_nodes.mod" "$(INTDIR)\local_smooth_fix.mod"\
 "$(INTDIR)\data_preparation.mod" "$(INTDIR)\update_curves.mod"\
 "$(INTDIR)\fix_reversed_curves.mod" "$(INTDIR)\form_discurves.mod"\
 "$(INTDIR)\disjoints_fix.mod"
   $(BuildCmds)

"$(INTDIR)\reset_discvs.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\update_left_n_right_states.mod"\
 "$(INTDIR)\produce_local_smooth_4_nodes.mod" "$(INTDIR)\local_smooth_fix.mod"\
 "$(INTDIR)\data_preparation.mod" "$(INTDIR)\update_curves.mod"\
 "$(INTDIR)\fix_reversed_curves.mod" "$(INTDIR)\form_discurves.mod"\
 "$(INTDIR)\disjoints_fix.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_RESET=\
	".\Debug\update_left_n_right_states.mod"\
	".\Debug\produce_local_smooth_4_nodes.mod"\
	".\Debug\local_smooth_fix.mod"\
	".\Debug\data_preparation.mod"\
	".\Debug\update_curves.mod"\
	".\Debug\fix_reversed_curves.mod"\
	".\Debug\form_discurves.mod"\
	".\Debug\disjoints_fix.mod"\
	
F90_MODOUT=\
	"reset_discvs"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\reset.obj" : $(SOURCE) $(DEP_F90_RESET) "$(INTDIR)"\
 "$(INTDIR)\update_left_n_right_states.mod"\
 "$(INTDIR)\produce_local_smooth_4_nodes.mod" "$(INTDIR)\local_smooth_fix.mod"\
 "$(INTDIR)\data_preparation.mod" "$(INTDIR)\update_curves.mod"\
 "$(INTDIR)\fix_reversed_curves.mod" "$(INTDIR)\form_discurves.mod"\
 "$(INTDIR)\disjoints_fix.mod"
   $(BuildCmds)

"$(INTDIR)\reset_discvs.mod" : $(SOURCE) $(DEP_F90_RESET) "$(INTDIR)"\
 "$(INTDIR)\update_left_n_right_states.mod"\
 "$(INTDIR)\produce_local_smooth_4_nodes.mod" "$(INTDIR)\local_smooth_fix.mod"\
 "$(INTDIR)\data_preparation.mod" "$(INTDIR)\update_curves.mod"\
 "$(INTDIR)\fix_reversed_curves.mod" "$(INTDIR)\form_discurves.mod"\
 "$(INTDIR)\disjoints_fix.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\interact\dis_jntn.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_DIS_J=\
	".\Release\solution.mod"\
	
F90_MODOUT=\
	"disjoints_fix"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\dis_jntn.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\disjoints_fix.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_DIS_J=\
	".\Debug/solution.mod"\
	
F90_MODOUT=\
	"disjoints_fix"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\dis_jntn.obj" : $(SOURCE) $(DEP_F90_DIS_J) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\disjoints_fix.mod" : $(SOURCE) $(DEP_F90_DIS_J) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\scan\sc_conls.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SC_CO=\
	".\Release\discontinuity_curves.mod"\
	
F90_MODOUT=\
	"scan_conservation_on_list"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\sc_conls.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\scan_conservation_on_list.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SC_CO=\
	".\Debug\discontinuity_curves.mod"\
	
F90_MODOUT=\
	"scan_conservation_on_list"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\sc_conls.obj" : $(SOURCE) $(DEP_F90_SC_CO) "$(INTDIR)"\
 "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

"$(INTDIR)\scan_conservation_on_list.mod" : $(SOURCE) $(DEP_F90_SC_CO)\
 "$(INTDIR)" "$(INTDIR)\discontinuity_curves.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\phys_flx\phy_updt\eva_sts.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_EVA_S=\
	".\Release\produce_polynomial_4_xx_yy.mod"\
	".\Release\euler_functions.mod"\
	
F90_MODOUT=\
	"evaluate_updated_xxyy"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\eva_sts.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\produce_polynomial_4_xx_yy.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

"$(INTDIR)\evaluate_updated_xxyy.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\produce_polynomial_4_xx_yy.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_EVA_S=\
	".\Debug\produce_polynomial_4_xx_yy.mod"\
	".\Debug\euler_functions.mod"\
	
F90_MODOUT=\
	"evaluate_updated_xxyy"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\eva_sts.obj" : $(SOURCE) $(DEP_F90_EVA_S) "$(INTDIR)"\
 "$(INTDIR)\produce_polynomial_4_xx_yy.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

"$(INTDIR)\evaluate_updated_xxyy.mod" : $(SOURCE) $(DEP_F90_EVA_S) "$(INTDIR)"\
 "$(INTDIR)\produce_polynomial_4_xx_yy.mod" "$(INTDIR)\euler_functions.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_oprt\prc_auxi\isl_fix.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_ISL_F=\
	".\Release\auxiliary_discontinuity_curves.mod"\
	".\Release\auxi_cc_prod.mod"\
	
F90_MODOUT=\
	"auxiliary_list_fix"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\isl_fix.obj" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_list_fix.mod" : $(SOURCE) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_ISL_F=\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\auxi_cc_prod.mod"\
	
F90_MODOUT=\
	"auxiliary_list_fix"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\isl_fix.obj" : $(SOURCE) $(DEP_F90_ISL_F) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"
   $(BuildCmds)

"$(INTDIR)\auxiliary_list_fix.mod" : $(SOURCE) $(DEP_F90_ISL_F) "$(INTDIR)"\
 "$(INTDIR)\auxiliary_discontinuity_curves.mod" "$(INTDIR)\auxi_cc_prod.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\interface\sh_sts_el_HR.f90

!IF  "$(CFG)" == "main_euler - Win32 Release"

NODEP_F90_SH_STS=\
	".\Release\solution.mod"\
	
F90_MODOUT=\
	"output_show_HR"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\sh_sts_el_HR.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\output_show_HR.mod" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "main_euler - Win32 Debug"

DEP_F90_SH_STS=\
	".\Debug/solution.mod"\
	
F90_MODOUT=\
	"output_show_HR"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\sh_sts_el_HR.obj" : $(SOURCE) $(DEP_F90_SH_STS) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

"$(INTDIR)\output_show_HR.mod" : $(SOURCE) $(DEP_F90_SH_STS) "$(INTDIR)"\
 "$(INTDIR)\solution.mod"
   $(BuildCmds)

!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
