# Abstract:
#	Template makefile for building Real-Time Windows Target compatible
#       real-time version of Simulink model using generated C code and
#       the built-in Open Watcom C/C++ Compiler.
#
# 	This makefile is designed to be used with GNU Make (gmake) which is
#       located in matlabroot/bin/win32.
#
# 	Note that this template is automatically customized by the Real-Time
#	Workshop build procedure to create "<model>.mk"
#
#       The following defines (macro names) can be used to modify the behavior
#       of the build:
#	  OPT_OPTS       - Optimization options.
#	  OPTS           - User-specified compiler options.
#         CPP_OPTS       - User-specified C++ compiler options.
#	  USER_SRCS      - Additional user objects, such as files needed by
#			   S-functions.
#	  USER_INCLUDES  - Additional include paths (i.e.
#			   "USER_INCLUDES=include-path1;include-path2")
#                          Use a '/' as a file separator instead of '\'.
#
#       Consider using the "Build process" dialog in Code Generation
#       options page instead of trying to change OPT_OPTS here.
#
#       This template makefile is designed to be used with a system target
#       file that contains 'rtwgensettings.ProjectDirSuffix', see grt.tlc .
#
# !!! THIS FILE IS AUTO-GENERATED !!!
# Copyright 1994-2011 The MathWorks, Inc.



#------------------------ Macros read by make_rtw -----------------------------
#
# The following macros are read by the code generation build procedure:
#  MAKECMD         - This is the command used to invoke the make utility
#  HOST            - What platform this template makefile is targeted for
#                    (i.e. PC or UNIX)
#  BUILD           - Invoke make from the code generation build procedure
#                    (yes/no)?
#  SYS_TARGET_FILE - Name of system target file.
#

MAKECMD         = "C:/PROGRA~1/MATLAB/R2012b/bin/win32/gmake"
HOST            = PC
BUILD           = yes
SYS_TARGET_FILE = rtwin.tlc
COMPILER_TOOL_CHAIN = default
MAKEFILE_FILESEP = /


#---------------------- Tokens expanded by make_rtw ---------------------------
#
# The following tokens, when wrapped with "|>" and "<|" are expanded by the
# code generation build procedure.
#
#  MODEL_NAME          - Name of the Simulink block diagram
#  MODEL_MODULES       - Any additional generated source modules
#  MAKEFILE_NAME       - Name of makefile created from template makefile <model>.mk
#  MATLAB_ROOT         - Path to were MATLAB is installed.
#  MATLAB_BIN          - Path to MATLAB executable.
#  S_FUNCTIONS         - List of S-functions.
#  S_FUNCTIONS_LIB     - List of S-functions libraries to link.
#  SOLVER              - Solver source file name
#  NUMST               - Number of sample times
#  TID01EQ             - yes (1) or no (0): Are sampling rates of continuous
#                        task (tid=0) and 1st discrete task equal.
#  NCSTATES            - Number of continuous states
#  BUILDARGS           - Options passed in at the command line.
#  MULTITASKING        - yes (1) or no (0): Is solver mode multitasking
#  EXT_MODE            - yes (1) or no (0): Build for external mode
#  EXTMODE_TRANSPORT   - Name of transport mechanism (e.g. tcpip, serial) for extmode
#  EXTMODE_STATIC      - yes (1) or no (0): Use static instead of dynamic mem alloc.
#  EXTMODE_STATIC_SIZE - Size of static memory allocation buffer.
#
#  CC_LISTING          - yes (1) or no (0): Generate assembly listings
#  REBUILD_ALL         - yes (1) or no (0): Rebuild all files

MODEL                := Control_linear
MODULES              := Control_linear_data.c rtGetInf.c rtGetNaN.c rt_nonfinite.c 
MAKEFILE             := Control_linear.mk
MATLAB_ROOT          := C:/PROGRA~1/MATLAB/R2012b
MATLAB_BIN           := C:/PROGRA~1/MATLAB/R2012b/bin
S_FUNCTIONS          := 
S_FUNCTIONS_LIB      := 
SOLVER               := 
NUMST                := 2
TID01EQ              := 1
NCSTATES             := 10
BUILDARGS            :=  GENERATE_REPORT=0 EXTMODE_STATIC_ALLOC=0 EXTMODE_STATIC_ALLOC_SIZE=1000000 TMW_EXTMODE_TESTING=0
MULTITASKING         := 0
EXT_MODE             := 1
EXTMODE_TRANSPORT    := 0

MODELREFS            := 
SHARED_SRC           := 
SHARED_SRC_DIR       := 
SHARED_BIN_DIR       := 
SHARED_LIB           := 
TARGET_LANG_EXT      := c
OPTIMIZATION_FLAGS   := 
ADDITIONAL_LDFLAGS   := 

# Real-Time Windows Target specific parameters
RTWINDIR             := C:/PROGRA~1/MATLAB/R2012b/toolbox/rtw/targets/rtwin
TARGETARCH           := win32
CC_LISTING           := 0
REBUILD_ALL          := 0

# ensure MATLAB_ROOT uses forward slashes - necessary for UNC paths
MATLAB_ROOT := $(subst \\,/,$(MATLAB_ROOT))

# some makefile magic
COMMA := ,
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)

# compute RXEXT based on target architecture
RXEXT := $(subst win,rxw,$(TARGETARCH))

#--------------------------- Model and reference models -----------------------
#
MODELLIB                  := Control_linearlib_rtwin.lib
MODELREF_LINK_LIBS        := 
MODELREF_LINK_RSPFILE     := Control_linear_ref.rsp
MODELREF_INC_PATH         := 
RELATIVE_PATH_TO_ANCHOR   := ..
MODELREF_TARGET_TYPE      := NONE


#--------------------------------- Tool Locations -----------------------------
#
SHELL   := cmd


#------------------------ External mode ---------------------------------------
#
# To add a new transport layer, see the comments in
#   <matlabroot>/toolbox/simulink/simulink/extmode_transports.m
ifeq ($(EXT_MODE),1)
  EXT_CC_OPTS := -DEXT_MODE
endif


#------------------------------ Include Path -----------------------------
#
# MATLAB includes
REQ_INCLUDES := $(MATLAB_ROOT)/simulink/include;$(MATLAB_ROOT)/extern/include;$(MATLAB_ROOT)/rtw/c/src;$(MATLAB_ROOT)/rtw/c/src/ext_mode/common

# target-specific and compiler-specific includes
REQ_INCLUDES := $(REQ_INCLUDES);$(RTWINDIR)/src

# additional includes
REQ_INCLUDES := $(REQ_INCLUDES);C:/Users/yk.L/Desktop/lsy/线性化PID/三通道（偏航串级期望至俯仰）/Control_linear_rtwin_win32;C:/Users/yk.L/Desktop/lsy/线性化PID/三通道（偏航串级期望至俯仰）

# shared includes
ifneq ($(SHARED_SRC_DIR),)
  REQ_INCLUDES := $(REQ_INCLUDES);$(SHARED_SRC_DIR)
endif

INCLUDES := $(USER_INCLUDES);.;$(RELATIVE_PATH_TO_ANCHOR);$(REQ_INCLUDES)$(MODELREF_INC_PATH)


# include platform-specific defines now when all template arguments are expanded
ifneq ($(wildcard $(RTWINDIR)/rtwin/rtwin_$(TARGETARCH).mk),)
  include $(RTWINDIR)/rtwin/rtwin_$(TARGETARCH).mk
else
  $(error Building of external mode executables is not yet supported on $(TARGETARCH))
endif


#-------------------------------- C Flags --------------------------------
#
OPT_OPTS ?= $(DEFAULT_OPT_OPTS)
COMMON_OPTS := $(REQ_OPTS) $(OPT_OPTS) $(OPTS) $(EXT_CC_OPTS)

REQ_DEFINES := -DUSE_RTMODEL -DMODEL=$(MODEL) -DRT -DNUMST=$(NUMST) \
               -DTID01EQ=$(TID01EQ) -DNCSTATES=$(NCSTATES) \
               -DMT=$(MULTITASKING) $(CC_REQ_DEFINES)

USER_INCLUDES ?=

# form the complete compiler command
CC += $(CC_REQ_OPTS) $(COMMON_OPTS) $(REQ_DEFINES)
CPP += $(CPP_REQ_OPTS) $(CPP_OPTS) $(COMMON_OPTS) $(REQ_DEFINES)


#------------------------------- Source Files ---------------------------------
#
# standard target
ifeq ($(MODELREF_TARGET_TYPE),NONE)
  PRODUCT     := $(RELATIVE_PATH_TO_ANCHOR)/$(MODEL).$(RXEXT)
  REQ_SRCS    := $(MODEL).$(TARGET_LANG_EXT) $(MODULES) \
                 rtwin_main.c rt_sim.c

ifeq ($(EXT_MODE),1)
  REQ_SRCS    += ext_svr.c updown_rtwin.c 
endif

# model reference target
else
  PRODUCT  := $(MODELLIB)
  REQ_SRCS := $(MODULES)
endif

SRCS := $(REQ_SRCS) $(USER_SRCS) $(S_FUNCTIONS)
SRCS += $(SOLVER)
OBJS := $(patsubst %.c,%.obj,$(patsubst %.cpp,%.obj,$(SRCS)))

SHARED_OBJS := $(addsuffix .obj, $(basename $(wildcard $(SHARED_SRC))))


#---------------------------- Additional Libraries ----------------------------
#
LIBS := 




#-------------------------- Rules ---------------------------------------
#
# decide what should get built
.DEFAULT_GOAL := $(if $(filter 1,$(REBUILD_ALL)), rebuild, $(PRODUCT))

# rule to rebuild everything unconditionally
.PHONY : rebuild
rebuild :
	$(MAKE) -B -f $(MAKEFILE) REBUILD_ALL=0

# libraries to link with the executable
ALLLIBS := $(strip $(LIBS) $(SHARED_LIB) $(MODELREF_LINK_LIBS))

# rules for linking the executable or modelref static library
ifeq ($(MODELREF_TARGET_TYPE),NONE)
$(PRODUCT) : $(OBJS) $(ALLLIBS)
	$(call LINKEXE, $@, $(OBJS), $(ALLLIBS))
	$(info ### Created Real-Time Windows Target module $(notdir $@))
else
$(PRODUCT) : $(OBJS) $(SHARED_LIB)
	$(call LIBEXE, $@, $(OBJS))
	$(info ### Created static library $@)
endif

# object build macros
CC_CPP := $(if $(filter cpp,$(TARGET_LANG_EXT)),$(CPP),$(CC))
define BUILDOBJ
	$(info ### Compiling $<)
	$(1) $(OUTFILE_OPT)"$@" "$<"
	$(CCLISTING)
endef

# rules for compiling objects
rtwin_main.obj : $(RTWINDIR)/src/rtwin_main.c $(MAKEFILE)
	$(call BUILDOBJ, $(CC_CPP))

%.obj : $(RTWINDIR)/src/%.c
	$(call BUILDOBJ, $(CC))

%.obj : $(MATLAB_ROOT)/rtw/c/src/ext_mode/common/%.c
	$(call BUILDOBJ, $(CC))

%.obj : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(call BUILDOBJ, $(CC))

%.obj : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	$(call BUILDOBJ, $(CPP))

%.obj : %.c
	$(call BUILDOBJ, $(CC))

%.obj : %.cpp
	$(call BUILDOBJ, $(CPP))

# additional sources
%.obj : $(MATLAB_ROOT)/rtw/c/src/%.c $(MAKEFILE)
	$(call BUILDOBJ, $(CC))

%.obj : $(MATLAB_ROOT)/rtw/c/src/ext_mode/common/%.c $(MAKEFILE)
	$(call BUILDOBJ, $(CC))

%.obj : $(MATLAB_ROOT)/rtw/c/src/rtiostream/rtiostreamtcpip/%.c $(MAKEFILE)
	$(call BUILDOBJ, $(CC))

%.obj : F:/Work/workforPenson/googoltech/GHP2002/GT400驱动程序/%.c $(MAKEFILE)
	$(call BUILDOBJ, $(CC))



%.obj : $(MATLAB_ROOT)/rtw/c/src/%.cpp $(MAKEFILE)
	$(call BUILDOBJ, $(CPP))

%.obj : $(MATLAB_ROOT)/rtw/c/src/ext_mode/common/%.cpp $(MAKEFILE)
	$(call BUILDOBJ, $(CPP))

%.obj : $(MATLAB_ROOT)/rtw/c/src/rtiostream/rtiostreamtcpip/%.cpp $(MAKEFILE)
	$(call BUILDOBJ, $(CPP))

%.obj : F:/Work/workforPenson/googoltech/GHP2002/GT400驱动程序/%.cpp $(MAKEFILE)
	$(call BUILDOBJ, $(CPP))



# simulink/src helper files
%.obj : $(MATLAB_ROOT)/simulink/src/%.c
	$(call BUILDOBJ, $(CC))

%.obj : $(MATLAB_ROOT)/simulink/src/%.cpp
	$(call BUILDOBJ, $(CPP))

# model-referencing shared objects
$(SHARED_BIN_DIR)/%.obj : $(SHARED_SRC_DIR)/%.c 
	$(call BUILDOBJ, $(CC))

$(SHARED_BIN_DIR)/%.obj : $(SHARED_SRC_DIR)/%.cpp 
	$(call BUILDOBJ, $(CPP))


# model-referencing shared library
$(SHARED_LIB) : $(SHARED_OBJS)
	$(call LIBEXE, $@, $^)
	$(info ### Created shared library $@)


# rules for building libraries



# rules for precompiled libraries

