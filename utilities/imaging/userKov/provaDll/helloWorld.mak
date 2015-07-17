# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

!IF "$(CFG)" == ""
CFG=helloWorld - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to helloWorld - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "helloWorld - Win32 Release" && "$(CFG)" !=\
 "helloWorld - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "helloWorld.mak" CFG="helloWorld - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "helloWorld - Win32 Release" (based on\
 "Win32 (x86) Dynamic-Link Library")
!MESSAGE "helloWorld - Win32 Debug" (based on\
 "Win32 (x86) Dynamic-Link Library")
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
# PROP Target_Last_Scanned "helloWorld - Win32 Debug"
RSC=rc.exe
MTL=mktyplib.exe
F90=fl32.exe

!IF  "$(CFG)" == "helloWorld - Win32 Release"

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

ALL : "$(OUTDIR)\helloWorld.dll"

CLEAN : 
	-@erase ".\Release\helloWorld.dll"
	-@erase ".\Release\prova2.obj"
	-@erase ".\Release\helloWorld.lib"
	-@erase ".\Release\helloWorld.exp"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo /MT
# ADD F90 /Ox /I "Release/" /c /nologo /MT
F90_PROJ=/Ox /I "Release/" /c /nologo /MT /Fo"Release/" 
F90_OBJS=.\Release/
# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE RSC /l 0x410 /d "NDEBUG"
# ADD RSC /l 0x410 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/helloWorld.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:windows /dll /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:windows /dll /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:windows /dll /incremental:no\
 /pdb:"$(OUTDIR)/helloWorld.pdb" /machine:I386 /out:"$(OUTDIR)/helloWorld.dll"\
 /implib:"$(OUTDIR)/helloWorld.lib" 
LINK32_OBJS= \
	"$(INTDIR)/prova2.obj"

"$(OUTDIR)\helloWorld.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "helloWorld - Win32 Debug"

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

ALL : "$(OUTDIR)\helloWorld.dll"

CLEAN : 
	-@erase ".\Debug\helloWorld.dll"
	-@erase ".\Debug\prova2.obj"
	-@erase ".\Debug\helloWorld.ilk"
	-@erase ".\Debug\helloWorld.lib"
	-@erase ".\Debug\helloWorld.exp"
	-@erase ".\Debug\helloWorld.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo /MT
# ADD F90 /Zi /I "Debug/" /c /nologo /MT
F90_PROJ=/Zi /I "Debug/" /c /nologo /MT /Fo"Debug/" /Fd"Debug/helloWorld.pdb" 
F90_OBJS=.\Debug/
# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE RSC /l 0x410 /d "_DEBUG"
# ADD RSC /l 0x410 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/helloWorld.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:windows /dll /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:windows /dll /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:windows /dll /incremental:yes\
 /pdb:"$(OUTDIR)/helloWorld.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/helloWorld.dll" /implib:"$(OUTDIR)/helloWorld.lib" 
LINK32_OBJS= \
	"$(INTDIR)/prova2.obj"

"$(OUTDIR)\helloWorld.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
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

# Name "helloWorld - Win32 Release"
# Name "helloWorld - Win32 Debug"

!IF  "$(CFG)" == "helloWorld - Win32 Release"

!ELSEIF  "$(CFG)" == "helloWorld - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\prova2.f90

"$(INTDIR)\prova2.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
