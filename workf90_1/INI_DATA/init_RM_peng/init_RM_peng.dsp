# Microsoft Developer Studio Project File - Name="init_RM" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=init_RM - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "init_RM_peng.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "init_RM_peng.mak" CFG="init_RM - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "init_RM - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "init_RM - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "init_RM - Win32 Release"

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
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "init_RM - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "init_RM - Win32 Release"
# Name "init_RM - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\adss_cell.f90
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\auxi_cuv\auxi_cc.f90
DEP_F90_AUXI_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\computational_information.mod"\
	".\Debug\critical_cells.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\AUXI_CUV\Auxi_cv1.f90
DEP_F90_AUXI_C=\
	".\Debug\auxiliary_critical_cells.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\reg_cuv\crit_cell.f90
DEP_F90_CRIT_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\crit_info.f90
DEP_F90_CRIT_I=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\solu\data_ext.f90
DEP_F90_DATA_=\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\reg_cuv\discontinuity_curve.f90
DEP_F90_DISCO=\
	".\Debug\critical_cells.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\main_itf\main_euler\euler\euler_functions.f90
# End Source File
# Begin Source File

SOURCE=..\..\MAIN_ITF\main_euler\EULER\eulernew.f90
DEP_F90_EULER=\
	".\Debug\data_extrapolation.mod"\
	".\Debug\euler_functions.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\auxi_cuv\fif_cell.f90
DEP_F90_FIF_C=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FORM_SL_MD.F90
DEP_F90_FORM_=\
	".\Debug\euler_functions.mod"\
	".\Debug\grid_map.mod"\
	".\Debug\numerical_integral_in_2D.mod"\
	".\Debug\out_intermediate.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\REST_CV\Fun_use.f90
DEP_F90_FUN_U=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\gif_cell.f90
DEP_F90_GIF_C=\
	".\Debug\grid.mod"\
	".\Debug\tools.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\grid.f90
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\grid_map.f90
DEP_F90_GRID_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\grid.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\init_RM_peng_modified.f90
DEP_F90_INIT_=\
	".\Debug\form_solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\node\nd_cls.f90
DEP_F90_ND_CL=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\node_plug_ring.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\node\nd_ring.f90
DEP_F90_ND_RI=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\Num_int.f90
DEP_F90_NUM_I=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\Num_int2.f90
DEP_F90_NUM_IN=\
	".\Debug\function_in_use.mod"\
	".\Debug\numerical_integrals.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\MAIN_ITF\main_euler\interface\outputs_euler1.f90
DEP_F90_OUTPU=\
	".\Debug\grid_map.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\dis_cuv\pif_cell.f90
DEP_F90_PIF_C=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\Ploynml.f90
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\sl_smth\sh_sths.f90
DEP_F90_SH_ST=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\solu\solution.f90
DEP_F90_SOLUT=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\grid_map.mod"\
	".\Debug\node_cells.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_strct\sl_smth\solution_smooth.f90
DEP_F90_SOLUTI=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	".\Debug\show_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\bs_oprt\tools.f90
DEP_F90_TOOLS=\
	".\Debug\polynomials.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
