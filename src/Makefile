F90 := ftn

EXE := gemx
ROOT := ..
BINDIR := $(ROOT)/bin
JOBDIR := $(ROOT)/jobs
SCRIPTDIR := $(ROOT)/scripts

SRCS := $(wildcard *.F90)
OBJS := $(SRCS:.F90=.o)
LIBS := -L${PETSC_PATH}/lib -lpetsc
INCS := -I${PETSC_PATH}/include
SHELLS := $(wildcard *.sh)
SCRIPTS := $(wildcard $(SCRIPTDIR)/*)

#Setup flags for various situations.
DEBUG ?= 0 #Use 'make DEBUG=1' to build with debug capabilities.
CPU ?= 0   #Use 'make CPU=1' to build for cpu.
OPT_COMMON := -r8 -Mfree $(INCS)
OPT_COMMON_REL := -O4
OPT_COMMON_DBG := -O0 -g -C -traceback -Mchkptr
OPT_CPU := -acc=host #Default to serial for cpu. -acc=multicore for parallel but should test against OpenMP on system.
OPT_GPU := -acc -Minfo=accel
OPT_CPU_DBG := -Mbounds #Note, bounds checking disabled on gpu by OpenACC.

ifeq ($(DEBUG),1)
	ifeq ($(CPU),1)
		JOB := job_dbg_cpu
		OPT :=  $(OPT_COMMON_DBG) $(OPT_COMMON) $(OPT_CPU) $(OPT_CPU_DBG)
	else
		JOB := job_dbg
		OPT :=  $(OPT_COMMON_DBG) $(OPT_COMMON) $(OPT_GPU)
	endif
else
	ifeq ($(CPU),1)
		JOB := job_cpu
		OPT := $(OPT_COMMON_REL) $(OPT_COMMON) $(OPT_CPU)
	else
		JOB := job
		OPT := $(OPT_COMMON_REL) $(OPT_COMMON) $(OPT_GPU)
	endif
endif

#If not on perlmutter use local run script and local compiler.
ifneq ($(NERSC_HOST),perlmutter)
	F90 := mpif90
	JOB := job_local
	SHELLS := $(filter-out env.sh,$(SHELLS))
endif

$(EXE): $(OBJS)
	$(F90) -o $(EXE) $(OPT) $(LIBS) $(OBJS)
	$(MAKE) install

pputil.o: pputil.F90
	$(F90) -c $(OPT) pputil.F90

gemx_com.o: gemx_com.F90 pputil.o
	$(F90) -c $(OPT) gemx_com.F90

equil.o: equil.F90 pputil.o gemx_com.o
	$(F90) -c $(OPT) equil.F90

gemx.o: gemx.F90 pputil.o gemx_com.o equil.o
	$(F90) -c $(OPT) gemx.F90

ionPush.o: ionPush.F90 gemx_com.o equil.o
	$(F90) -c $(OPT) ionPush.F90

outd.o: outd.F90 gemx_com.o equil.o
	$(F90) -c $(OPT) outd.F90

fcnt.o: fcnt.F90
	$(F90) -c $(OPT) fcnt.F90

clean:
	$(RM) *.o *.lst *.mod gemx
	$(RM) -r $(BINDIR)

install:
	mkdir -p $(BINDIR) $(BINDIR)/out $(BINDIR)/dump $(BINDIR)/matrix
	cp -f gemx $(BINDIR)
	cp -f gemx.in $(BINDIR)
	cp -f *.dat $(BINDIR)
	cp -f $(JOBDIR)/$(JOB) $(BINDIR)
	cp -f $(SCRIPTS) $(BINDIR)
	cp -f $(SHELLS) $(BINDIR)