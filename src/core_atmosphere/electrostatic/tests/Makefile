FC ?= mpifort
FFLAGS ?= -O0 -g -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8

MPAS_ROOT := ../../../..
INCLUDES := -I$(MPAS_ROOT)/src/framework -I$(MPAS_ROOT)/src/operators -I..
LIBS := $(MPAS_ROOT)/src/operators/libops.a $(MPAS_ROOT)/src/framework/libframework.a
ELECTROSTATIC_OBJS = ../mpas_electrostatic_bcs.o
MMS_OBJS = ../mpas_electrostatic_source.o ../mpas_electrostatic_benchmarks.o

TESTS = test_elliptic_operator test_pcg_solver test_mms_source test_e_reconstruction

all: $(TESTS)

test_elliptic_operator: test_elliptic_operator.F90 $(LIBS)
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ test_elliptic_operator.F90 $(LIBS)

test_pcg_solver: test_pcg_solver.F90 $(LIBS)
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ test_pcg_solver.F90 $(LIBS)

test_mms_source: test_mms_source.F90
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ test_mms_source.F90 $(MMS_OBJS) $(LIBS)

test_e_reconstruction: test_e_reconstruction.F90 $(ELECTROSTATIC_OBJS) $(LIBS)
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ test_e_reconstruction.F90 $(ELECTROSTATIC_OBJS) $(LIBS)

clean:
	$(RM) $(TESTS) *.o *.mod
