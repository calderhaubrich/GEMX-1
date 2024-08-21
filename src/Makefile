F90  = mpif90

SRCS = gem_com.F90 equil.F90 gem.F90 outd.F90 fcnt.F90
OBJS = gem_com.o equil.o gem.o ionPush.o outd.o fcnt.o
PLIB = pputil.o
LIBS = -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc
INCS = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

DEBUG ?= 1 #Default to off. Use 'make DEBUG=1' to build in debug.
OPT_COMMON = -acc -r8 -Mbounds -Mfree -Minfo=accel $(INCS)
ifeq ($(DEBUG),1)
	OPT = -O0 -g -C -traceback -Mchkptr $(OPT_COMMON)
else
	OPT = -O4 $(OPT_COMMON)
endif

gemx: equil.o gem.o ionPush.o outd.o fcnt.o pputil.o gem_com.o
	$(F90) -o gemx $(OPT) $(LIBS) $(OBJS) $(PLIB)

pputil.o: pputil.F90
	$(F90) -c $(OPT) pputil.F90

gem_com.o: gem_com.F90 pputil.o
	$(F90) -c $(OPT) gem_com.F90

equil.o: equil.F90 pputil.o gem_com.o
	$(F90) -c $(OPT)  equil.F90

gem.o: gem.F90 pputil.o gem_com.o equil.o
	$(F90) -c $(OPT) gem.F90

ionPush.o: ionPush.F90 gem_com.o equil.o
	$(F90) -c $(OPT) ionPush.F90

outd.o: outd.F90 gem_com.o equil.o
	$(F90) -c $(OPT) outd.F90

fcnt.o: fcnt.F90
	$(F90) -c $(OPT) fcnt.F90

clean:
	rm -f *.o *.lst *.mod gemx out/* dump/* xpp matrix/* testne testphi psi_test run.out flag_debug debug.dat test* mask* gemx_timing.txt