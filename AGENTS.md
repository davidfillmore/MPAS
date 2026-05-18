# AGENTS.md

## Session Startup

- Read `CLAUDE.md` first. It points to `../MPAS-Papers/CLAUDE.md`, which holds the shared MPAS workflow, planning, and journal rules.
- Plans and session history live in `../MPAS-Papers/`, not in this public code fork.

## Build Discipline

- The `mpas` conda environment on this host is for Python helpers
  (`netCDF4`, `numpy`, plotting, runner scripts). It does not provide the
  Fortran MPI wrapper used to build MPAS here.
- On the local macOS development host, use the LLVM build stack: `clang`, `flang`, and Homebrew Open MPI wrappers configured with:

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export OMPI_FC=flang
export OMPI_CC=clang
export OMPI_CXX=clang++
```

- Do not mix compiler families. `flang` `.mod` files are incompatible with `gfortran` `.mod` files. If the tree was built with `llvm`, standalone Fortran tests must also be compiled through `mpifort` with `OMPI_FC=flang`.
- For the atmosphere core on this host, use:

```bash
make -j4 llvm CORE=atmosphere PIO="$HOME/software" NETCDF=/opt/homebrew PNETCDF="$HOME/software" PRECISION=double
```

- For Poisson electrostatic standalone Fortran tests after an LLVM build, use:

```bash
OMPI_FC=flang OMPI_CC=clang OMPI_CXX=clang++ \
  make -C src/core_atmosphere/electrostatic mpas_electrostatic_source.o

OMPI_FC=flang OMPI_CC=clang OMPI_CXX=clang++ \
  make -C src/core_atmosphere/electrostatic/tests FC=mpifort test_mms_source
```

## Run Discipline

- Generated MPAS run artifacts, logs, NetCDF output, and plots belong under `~/Data/MPAS/`, not in this repository.
- On the local macOS PIO build, every MPAS stream used in multi-rank runs should set `io_type="netcdf"` explicitly.
- The synthetic tripole supercell runner is
  `src/core_atmosphere/electrostatic/scripts/run_tripole_supercell.py`; its
  default output root is `~/Data/MPAS/poisson_tripole_supercell/`.
