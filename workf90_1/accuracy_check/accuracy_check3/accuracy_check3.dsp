# Microsoft Developer Studio Project File - Name="accuracy_3" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=accuracy_3 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "accuracy_check3.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "accuracy_check3.mak" CFG="accuracy_3 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "accuracy_3 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "accuracy_3 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "accuracy_3 - Win32 Release"

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

!ELSEIF  "$(CFG)" == "accuracy_3 - Win32 Debug"

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

# Name "accuracy_3 - Win32 Release"
# Name "accuracy_3 - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\accuracy_check3.f90
NODEP_F90_ACCUR=\
	".\Debug\input_setup.mod"\
	".\Debug\input_setup_f.mod"\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\adss_cell.f90
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\AUXI_CUV\Auxi_cc.f90
DEP_F90_AUXI_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\computational_information.mod"\
	".\Debug\critical_cells.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\AUXI_CUV\Auxi_cv.f90
DEP_F90_AUXI_C=\
	".\Debug\auxiliary_critical_cells.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\REG_CUV\crit_cell.f90
NODEP_F90_CRIT_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\crit_info.f90
NODEP_F90_CRIT_I=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\geo_info_cell.mod"\
	".\Debug\phy_info_cell.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\SOLU\Data_ext.f90
DEP_F90_DATA_=\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\REG_CUV\discontinuity_curve.f90
NODEP_F90_DISCO=\
	".\Debug\critical_cells.mod"\
	".\Debug\grid.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\MAIN_ITF\main_euler\EULER\euler_functions.f90
# End Source File
# Begin Source File

SOURCE=..\..\MAIN_ITF\main_euler\EULER\Eulern.f90
NODEP_F90_EULER=\
	".\Debug\data_extrapolation.mod"\
	".\Debug\euler_functions.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\AUXI_CUV\Fif_cell.f90
DEP_F90_FIF_C=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\Gif_cell.f90
NODEP_F90_GIF_C=\
	".\Debug\grid.mod"\
	".\Debug\tools.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\Grid.f90
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\Grid_map.f90
NODEP_F90_GRID_=\
	".\Debug\adss_info_cell.mod"\
	".\Debug\auxiliary_discontinuity_curves.mod"\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\grid.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\inputf.f90
NODEP_F90_INPUT=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\MAIN_ITF\main_euler\interface\inputs_euler.f90
DEP_F90_INPUTS=\
	".\Debug\solution.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\NODE\Nd_cls.f90
NODEP_F90_ND_CL=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\node_plug_ring.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\NODE\Nd_ring.f90
NODEP_F90_ND_RI=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\DIS_CUV\Pif_cell.f90
NODEP_F90_PIF_C=\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\polynomial.f90
# End Source File
# Begin Source File

SOURCE=..\..\MAIN_ITF\main_euler\EULER\Sh_sths.f90
DEP_F90_SH_ST=\
	".\Debug\euler_functions.mod"\
	".\Debug\physical_state.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\SOLU\Solution.f90
NODEP_F90_SOLUT=\
	".\Debug\discontinuity_curves.mod"\
	".\Debug\grid_map.mod"\
	".\Debug\node_cells.mod"\
	".\Debug\solu_in_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_STRCT\SL_SMTH\solution_smooth.f90
NODEP_F90_SOLUTI=\
	".\Debug\grid.mod"\
	".\Debug\physical_state.mod"\
	".\Debug\show_smth.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\BS_OPRT\Tools.f90
NODEP_F90_TOOLS=\
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
