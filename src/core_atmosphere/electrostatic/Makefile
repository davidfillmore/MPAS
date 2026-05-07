.SUFFIXES: .F .o

OBJS = mpas_electrostatic_benchmarks.o \
       mpas_electrostatic_source.o

all: $(OBJS)

clean:
	$(RM) *.o *.mod *.f90
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators $(MPAS_ESMF_INC)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators $(MPAS_ESMF_INC)
endif
