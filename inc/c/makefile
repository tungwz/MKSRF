# FILE:
#     makefile
#
# PURPOSE:
#     build the avh2hdf program 

#**************************************************************
# Compiler
#**************************************************************
CC          = icc

#**************************************************************
# Directories
#**************************************************************
SRC_DIR     = . 

SRCDIR      += $(SRC_DIR)
VPATH       = $(SRCDIR)

#**************************************************************
# Compiler options
#**************************************************************

INCLUDE     = -I$(HDFINC) -I$(HDFEOS_INC) -I$(HDFEOS_HOME)/gctp/include 
LDFLAGS     = -L$(HDFEOS_LIB) -lhdfeos -lGctp -L$(HDFLIB) -lmfhdf -ldf -ljpeg -lz -lm -openmp
CC_OPT      = -g -c -O2 -openmp $(INCLUDE)

#**************************************************************
# The name of executable file
#**************************************************************
EXE_FILE    = test

#**************************************************************
# Files to be compiled
#**************************************************************
SRC_C       = $(foreach dir, $(SRCDIR), $(wildcard $(dir)/*.c))

#**************************************************************
# All link files for the object
#**************************************************************
OBJ_C       = $(patsubst %.c, %.o, $(SRC_C))
#OBJ_C       = $(notdir $(patsubst %.c, %.o, $(SRC_C)))

#**************************************************************
# Rules for make
#**************************************************************
$(EXE_FILE):$(OBJ_C)
	$(CC) -o $@ $^ $(LDFLAGS) 
	@echo $(EXE_FILE) "is now up2dated!"

$(OBJ_C):%.o:%.c
	$(CC) $(CC_OPT) -o $@ $<

clean:
	rm -f *.o
