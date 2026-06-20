.SUFFIXES: .F .o

# GNU make has a built-in FC=f77; use a normal assignment so local standalone
# object builds default to mpifort while command-line FC=... overrides still work.
FC = mpifort
FFLAGS ?= -O0 -g -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8

OBJS = mpas_electrostatic_benchmarks.o \
       mpas_electrostatic_bcs.o \
       mpas_electrostatic_driver.o \
       mpas_electrostatic_source.o \
       mpas_electrostatic_vertical_remap.o

all: $(OBJS)

clean:
	$(RM) *.o *.mod *.f90
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../framework -I../../operators $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../framework -I../../operators $(MPAS_ESMF_INC)
endif
