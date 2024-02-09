#
# Makefile for VANILLA version 1.0
#
#=================================================
# uncomment the following to enable OpenMP support
#OMP_FLAG =-fopenmp
#=================================================
#===== Macro definitions =====

SHELL = /bin/sh
VPATH = .:../include
VERSION = VANILLA_1.0
ENGINE = VANILLA10

EXE =
FFM_LIB =

buildat=$(shell LANG=C date)
#buildat=$(shell LANG=C date +'%a, %_d %b %Y %H:%M:%S %z')
BUILDDATA = -DBUILDAT='"$(buildat)"' -DBUILDON='"$(shell hostname)"'

#### default setting ##########

CXX = g++
LINKS	=  -lm -lstdc++ 
#LINKS	=  -lm -lfftw3
MAKE = make

#PF_ENGINEARCH=linux_64
PF_ENGINEARCH=win64

PF_FILES=/cygdrive/c/OCTA8.3/GOURMET
PF_ENGINE=/cygdrive/c/OCTA8.3/ENGINES

TARGET_DIR = $(PF_ENGINE)/bin/$(PF_ENGINEARCH)
#OBJ_DIR = objects
OBJ_DIR = objects/$(PF_ENGINEARCH)
PFLIB_DIR = $(PF_FILES)/lib/$(PF_ENGINEARCH)
PFENGINELIB_DIR = $(PF_ENGINE)/lib/$(PF_ENGINEARCH)

PFINC_DIR = $(PF_FILES)/include/
PFLIB = $(PFLIB_DIR)/libplatform.a


INCLUDE	= -I../include \
	   	-I$(PFINC_DIR) 

#### Platform dependent setting #######

# for Win32(MinGW)
ifeq ($(PF_ENGINEARCH),win32)
    CC = i686-w64-mingw32-gcc
    CXX = i686-w64-mingw32-g++
# The following options usually gives better performance
#	OPT_FLAG = -O3 -march=native -DMINGW
#	OPT_FLAG = -O3 -mtune=generic -DMINGW
	OPT_FLAG = -O3 -DMINGW
	LDFLAGS  = -static
	EXE = .exe
#	MIN_LIBS =
#	MIN_SRC = 
endif

# for Win64(MinGW)
ifeq ($(PF_ENGINEARCH),win64)
    CC = x86_64-w64-mingw32-gcc
    CXX = x86_64-w64-mingw32-g++
# The following options usually gives better performance
#	OPT_FLAG = -O3 -march=native -DMINGW
#	OPT_FLAG = -O3 -mtune=generic -DMINGW
	OPT_FLAG = -O3 -DMINGW
	LDFLAGS  = -static
	EXE = .exe
#	MIN_LIBS =
#	MIN_SRC = 
endif

# for Cygwin gcc
ifeq ($(PF_ENGINEARCH),cygwin)
# The following options usually gives better performance
#	OPT_FLAG = -O3 -march=native -DMINGW
#	OPT_FLAG = -O3 -mtune=generic -DMINGW
	OPT_FLAG = -O3 -DMINGW
	LDFLAGS  = -static
	EXE = .exe
##	MIN_LIBS =
##	MIN_SRC = 
endif

#for x86 linux
ifeq ($(PF_ENGINEARCH),linux)
# The following options usually gives better performance
#	OPT_FLAG = -O3 -march=native -DMINGW
#	OPT_FLAG = -O3 -mtune=generic -DMINGW
	OPT_FLAG = -O3 -DLINUX
#	LDFLAGS = -static
endif

#for x86_64 linux
ifeq ($(PF_ENGINEARCH),linux_64)
# This option usually gives better performance
#	OPT_FLAG = -O3 -march=native -DMINGW
#	OPT_FLAG = -O3 -mtune=generic -DMINGW
	OPT_FLAG = -O3 -DLINUX
#	LDFLAGS = -static
endif

#for Mac OS X
ifeq ($(PF_ENGINEARCH),macosx)
	OPT_FLAG =-O3
endif

#for AIX(SMP)
ifeq ($(PF_ENGINEARCH),aix)
	LINKS =
#	OPT_FLAG = -O4
#	OPT_FLAG = -O3 -qarch=auto -qtune=auto
	OPT_FLAG = -O3 -qsmp=omp:noauto -qstrict -A32767 -qarch=auto -qtune=auto -qtemplateregistry -DAIX
	LDFLAGS = -Wl,-blpdata -qsmp=omp
	CXX =xlC_r
## ... make must be GNU make !!
#	MAKE    = gmake
#	MIN_LIBS =
#	MIN_SRC =
endif

#########################################

MAKEDEPEND = /usr/bin/makedepend -Y

VANILLA_BASE_SRCS = \
		DataIO.cpp \
		Field.cpp \
		FreeEnergy.cpp \
		MyTools.cpp \
		Stress.cpp \
		Vector3c.cpp \
		VANILLA_main.cpp \
		Vector3d.cpp \



VANILLA_SRCS = $(VANILLA_BASE_SRCS)

OBJS_TMP = $(patsubst %.cpp,%.o,$(VANILLA_SRCS))
OBJS = $(addprefix $(OBJ_DIR)/,$(OBJS_TMP))

MISC_LIBS = $(MIN_LIBS) $(FFM_LIB)

APILIB_OPT = -DUSE_COMMON_SUPORT 
DSFMT_FLAG =  -DDSFMT_MEXP=19937

#ifeq ($(CXX),g++)
#	OTHER_FLAGS = -Wno-write-strings
#endif
#CXXFLAGS	= $(OPT_FLAG) $(OMP_FLAG) $(OTHER_FLAGS) $(APILIB_OPT) $(INCLUDE)
CXXFLAGS	= $(OPT_FLAG) $(OMP_FLAG) $(APILIB_OPT) $(INCLUDE) $(DSFMT_FLAG)

#<<<  Create a executable module. >>>
all :
	if [ ! -d $(OBJ_DIR) ]; then  \
		 mkdir -p $(OBJ_DIR); \
	fi
	$(MAKE) $(OBJ_DIR)/$(ENGINE)
$(OBJ_DIR)/$(ENGINE)	: $(OBJS) $(LIBS)
	$(CXX) -o $(OBJ_DIR)/$(ENGINE) $(LDFLAGS) $(OMP_FLAG) $(OBJS) \
	$(PFLIB) $(MISC_LIBS) $(LINKS)

#===== rules for VANILLA objs =====
$(OBJ_DIR)/%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 
$(OBJ_DIR)/VANILLA_main.o : VANILLA_main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ VANILLA_main.cpp


clean	:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/$(PF_ENGINEARCH)$(EXE) $(OBJ_DIR)/core

install	:
	if [ ! -d $(TARGET_DIR) ]; then  \
		 mkdir -p $(TARGET_DIR); \
	fi
	if [ -f $(OBJ_DIR)/$(ENGINE)$(EXE) ]; then \
		cp -f $(OBJ_DIR)/$(ENGINE)$(EXE) $(TARGET_DIR)/$(ENGINE)$(EXE);\
	fi

depend :
	$(MAKEDEPEND) -p$(OBJ_DIR)/ -- $(CXXFLAGS) -- $(VANILLA_SRCS) 2>/dev/null

#===== dependency created by makedepend (optional) =====
# DO NOT DELETE THIS LINE -- make depend depends on it.
