# FILE:
#     makefile
#
# PURPOSE:
#     build the tools program

#**************************************************************
# Compiler
#**************************************************************
CC          = icc
FC          = gfortran

#**************************************************************
# Directories
#**************************************************************
INC_DIR     = ./inc/c
INC_DIR     += ./inc/for
SRC_DIR     = .

SRCDIR      += $(INC_DIR) $(SRC_DIR)
VPATH       = $(SRCDIR)

#**************************************************************
# Compiler options
#**************************************************************

INCLUDE     = -I/opt/netcdf/include -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_HOME)/gctp/include -I./inc/c -I./inc/for -I/usr/include
LDFLAGS     = -L$(HDFEOS_LIB) -lhdfeos -lGctp -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lm -qopenmp -L/usr/lib -lnetcdf
CC_OPT      = -g -c -qopenmp $(INCLUDE)
FC_OPT      = -g -c -qopenmp $(INCLUDE)

#**************************************************************
# The name of executable file
#**************************************************************
EXE1       = makereglc
EXE2       = makereglai
EXE3       = makeregvcf
EXE4       = makeglobal15slc
EXE5       = makeregyearavglai

#**************************************************************
# Files to be compiled
#**************************************************************
SRC_C       = $(foreach dir, $(SRC_DIR), $(wildcard $(dir)/*.c))
SRC_F       = $(foreach dir, $(SRC_DIR), $(wildcard $(dir)/*.f90))
INC_C       = $(foreach dir, $(INC_DIR), $(wildcard $(dir)/*.c))
INC_F       = $(foreach dir, $(INC_DIR), $(wildcard $(dir)/*.f90))

#**************************************************************
# All link files for the object
#**************************************************************
OBJ_C       = $(patsubst %.c, %.o, $(SRC_C))
OBJ_F       = $(patsubst %.f90, %.o, $(SRC_F))
OBJ_INC_C   = $(patsubst %.c, %.o, $(INC_C))
OBJ_INC_F   = $(patsubst %.f90, %.o, $(INC_F))

#**************************************************************
# Rules for make
#**************************************************************
all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5)

$(EXE1):MakeReg_LC.o $(OBJ_INC_C)
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo $(EXE1) "is now up2dated!"
	@echo "!Warning: Before running, define DTYPE first in modis.h."
	@echo ""

$(EXE2):MakeReg_LAI.o $(OBJ_INC_C)
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo $(EXE2) "is now up2dated!"
	@echo "!Warning: Before running, define DTYPE first in modis.h."
	@echo ""

$(EXE3):MakeReg_VCF.o $(OBJ_INC_C)
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo $(EXE3) "is now up2dated!"
	@echo "!Warning: Before running, define DTYPE first in modis.h."
	@echo ""

$(EXE4):MakeGlobal15s_LC.o $(OBJ_INC_C)
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo $(EXE4) "is now up2dated!"
	@echo "!Warning: Before running, define DTYPE first in modis.h."
	@echo ""

$(EXE5):MakeRegYearAvg_LAI.o $(OBJ_INC_C)
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo ""
	@echo $(EXE5) "is now up2dated!"
	@echo "!Warning: Before running, define DTYPE first in modis.h."
	@echo ""


$(OBJ_C):%.o:%.c
	$(CC) $(CC_OPT) -o $@ $<

$(OBJ_F):%.o:%.f90
	$(FC) $(FC_OPT) -o $@ $<

$(OBJ_INC_C):%.o:%.c
	$(CC) $(CC_OPT) -o $@ $<

$(OBJ_INC_F):%.o:%.f90
	$(FC) $(FC_OPT) -o $@ $<

$(OBJ_C): tools.h
$(OBJ_INC_C): modis.h

clean:
	rm -f *.o *.mod
