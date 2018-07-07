# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=cycles - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to cycles - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "cycles - Win32 Release" && "$(CFG)" != "cycles - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "cycles.mak" CFG="cycles - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cycles - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "cycles - Win32 Debug" (based on "Win32 (x86) Console Application")
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
# PROP Target_Last_Scanned "cycles - Win32 Debug"
RSC=rc.exe
F90=fl32.exe

!IF  "$(CFG)" == "cycles - Win32 Release"

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

ALL : "$(OUTDIR)\cycles.exe"

CLEAN : 
	-@erase ".\Release\cycles.exe"
	-@erase ".\Release\cycl.obj"
	-@erase ".\Release\grid.mod"
	-@erase ".\Release\grid.obj"

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
BSC32_FLAGS=/nologo /o"$(OUTDIR)/cycles.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/cycles.pdb" /machine:I386 /out:"$(OUTDIR)/cycles.exe" 
LINK32_OBJS= \
	".\Release\cycl.obj" \
	".\Release\grid.obj"

"$(OUTDIR)\cycles.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "cycles - Win32 Debug"

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

ALL : "$(OUTDIR)\cycles.exe"

CLEAN : 
	-@erase ".\Debug\cycles.exe"
	-@erase ".\Debug\grid.obj"
	-@erase ".\Debug\cycl.obj"
	-@erase ".\Debug/grid.mod"
	-@erase ".\Debug\cycles.ilk"
	-@erase ".\Debug\cycles.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Zi /I "Debug/" /c /nologo
F90_PROJ=/Zi /I "Debug/" /c /nologo /Fo"Debug/" /Fd"Debug/cycles.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/cycles.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/cycles.pdb" /debug /machine:I386 /out:"$(OUTDIR)/cycles.exe" 
LINK32_OBJS= \
	".\Debug\grid.obj" \
	".\Debug\cycl.obj"

"$(OUTDIR)\cycles.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

# Name "cycles - Win32 Release"
# Name "cycles - Win32 Debug"

!IF  "$(CFG)" == "cycles - Win32 Release"

!ELSEIF  "$(CFG)" == "cycles - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=\workf90_1\bs_strct\grid.f90

!IF  "$(CFG)" == "cycles - Win32 Release"

F90_MODOUT=\
	"grid"


BuildCmds= \
	$(F90) $(F90_PROJ) $(SOURCE) \
	

"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

"$(INTDIR)\grid.mod" : $(SOURCE) "$(INTDIR)"
   $(BuildCmds)

!ELSEIF  "$(CFG)" == "cycles - Win32 Debug"

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

SOURCE=.\cycl.f90

!IF  "$(CFG)" == "cycles - Win32 Release"

NODEP_F90_CYCL_=\
	".\Release\grid.mod"\
	

"$(INTDIR)\cycl.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\grid.mod"


!ELSEIF  "$(CFG)" == "cycles - Win32 Debug"

DEP_F90_CYCL_=\
	".\Debug/grid.mod"\
	

"$(INTDIR)\cycl.obj" : $(SOURCE) $(DEP_F90_CYCL_) "$(INTDIR)"\
 "$(INTDIR)\grid.mod"


!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
