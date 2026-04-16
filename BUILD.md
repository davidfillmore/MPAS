# MPAS Build System

This document describes the build systems available for MPAS and their configuration options.

## Build Systems

MPAS supports two build systems:

1. **CMake** (recommended for new users)
2. **Legacy Makefile** (platform-specific targets)

## CMake Build

### Basic Usage

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `DO_PHYSICS` | ON | Enable built-in physics schemes |
| `MPAS_DOUBLE_PRECISION` | ON | Use 64-bit floating point |
| `MPAS_PROFILE` | OFF | Enable GPTL profiling |
| `MPAS_OPENMP` | OFF | Enable OpenMP parallelization |
| `BUILD_SHARED_LIBS` | ON | Build shared libraries |
| `MPAS_USE_PIO` | OFF | Use Parallel I/O library |

### Selecting Cores

By default, both `atmosphere` and `init_atmosphere` cores are built. To build specific cores:

```bash
cmake -DCORES="atmosphere" ..
```

### Example Configurations

**Basic atmosphere build:**
```bash
cmake -DCORES="atmosphere;init_atmosphere" \
      -DMPAS_DOUBLE_PRECISION=ON \
      -DDO_PHYSICS=ON \
      ..
```

**Debug build:**
```bash
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DMPAS_PROFILE=ON \
      ..
```

## Legacy Makefile Build

The `Makefile` supports various platform-specific targets:

### Available Targets

| Target | Compiler Suite |
|--------|---------------|
| `gnu` | GNU Fortran/C/C++ (traditional) |
| `gfortran` | GNU Fortran/C/C++ (conda-forge style) |
| `llvm` | LLVM flang/clang/clang++ |
| `xlf` | IBM XL compilers |
| `xlf-summit-omp-offload` | IBM XL with OpenMP offload (Summit) |
| `ftn` | Cray compilers |
| `intel` | Intel compilers |
| `pgi` | PGI/NVIDIA compilers |
| `nag` | NAG Fortran |

### Usage

```bash
make -j8 gfortran CORE=atmosphere \
  PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double

make -j8 gfortran CORE=init_atmosphere \
  PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```

**Note:** Use `-j8` (or appropriate number for your system) for parallel compilation.

### GCC Build on Ubuntu (conda-forge)

Building MPAS with GCC/gfortran on Ubuntu using a conda-forge toolchain.
This is the build path used on the development host in this repo — the
system `gcc`/`mpi` are **not** used. All compilers and libraries come
from the `mpas` conda environment except for PIO, which is built from
source and installed under `$HOME/software`.

#### Prerequisites

Create a conda environment with the full toolchain:
```bash
conda create -n mpas -c conda-forge \
  gcc gxx gfortran cmake openmpi libnetcdf netcdf-fortran libpnetcdf \
  libblas liblapack pkg-config make
conda activate mpas
```

Verify the environment resolves to the conda-forge toolchain (expected
values on this host):

| Check | Expected |
|-------|----------|
| `gfortran --version` | `GNU Fortran (conda-forge gcc 15.2.0-18) 15.2.0` |
| `mpifort --version` | wraps the same conda-forge gfortran |
| `which mpifort` | `$CONDA_PREFIX/bin/mpifort` |
| `nc-config --prefix` | `$CONDA_PREFIX` (i.e. `~/miniconda3/envs/mpas`) |
| `pkg-config --libs netcdf-fortran` | resolves without error |

If `which gfortran` points anywhere outside `$CONDA_PREFIX/bin`, the
conda env is not active or is shadowed by a `PATH` entry earlier in
`.bashrc` — fix that before building.

#### PIO Library

PIO is not available via conda and must be built from source. Install to
`$HOME/software` (the location the MPAS build expects in `PIO=`):

```bash
git clone https://github.com/NCAR/ParallelIO.git
cd ParallelIO && mkdir build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/software \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DPIO_ENABLE_TIMING=OFF
make -j8 && make install
```

Verify the install landed:
```bash
ls $HOME/software/lib/libpio*
# expected: libpioc.a  libpiof.a  libpio.settings
```

#### Building MPAS

From the MPAS repo root, with the `mpas` conda env active:

```bash
conda activate mpas
export NETCDF="$CONDA_PREFIX"
export PNETCDF="$CONDA_PREFIX"
export PIO="$HOME/software"

make -j8 gfortran \
  CORE=atmosphere \
  PIO="$PIO" \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double
```

For `init_atmosphere` (run as a second invocation — the Makefile builds
one core per invocation):
```bash
make -j8 gfortran \
  CORE=init_atmosphere \
  PIO="$PIO" \
  NETCDF="$NETCDF" \
  PNETCDF="$PNETCDF" \
  PRECISION=double
```

A successful build ends with the MPAS banner. On this host it prints:

```
MPAS was built with default double-precision reals.
Debugging is off.
Parallel version is on.
Using the mpi_f08 module.
...
MPAS was not linked with the MUSICA-Fortran library.
MPAS was NOT linked with the Scotch graph partitioning library.
...
Using the PIO 2.x library.
MPAS was built with the embedded ESMF timekeeping library.
```

`atmosphere_model` (≈9 MB) is left in the repo root.

#### One-liner rebuild (this host)

For convenience, the full command used on this host:

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate mpas && \
  export NETCDF="$CONDA_PREFIX" PNETCDF="$CONDA_PREFIX" PIO="$HOME/software" && \
  make -j8 gfortran CORE=atmosphere \
    PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
```

#### GCC Build Technical Details

The `gfortran` target configures:
- **MPI Wrappers:** Uses `mpif90`/`mpicc`/`mpicxx` (gfortran is the default wrapper compiler)
- **Fortran Flags:** `-fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-form -ffree-line-length-none`
- **MPI Module:** Uses `mpi_f08` module natively (no `NOMPIMOD` needed)

#### GCC-Specific Notes

1. **Conda environment isolation.** All builds (PIO, MPAS) must happen within
   the same conda environment to ensure consistent compiler/library paths.

2. **Physics data files.** The `checkout_data_files.sh` script in
   `src/core_atmosphere/physics/` downloads WRF lookup tables from
   [MPAS-Data](https://github.com/MPAS-Dev/MPAS-Data) on first build.

### LLVM Build on macOS (Homebrew)

Building MPAS with LLVM compilers (flang/clang) on macOS requires special configuration due to compiler and library compatibility issues.

#### Prerequisites

Install LLVM compilers and dependencies via Homebrew:
```bash
brew install llvm flang open-mpi netcdf netcdf-fortran
```

**Important:** Homebrew's Open MPI is compiled with gfortran, so its Fortran module files (`.mod`) are incompatible with flang. The build system handles this via the `NOMPIMOD` flag.

#### PnetCDF Requirement

Homebrew's PnetCDF is built with gfortran, making its Fortran interface incompatible with flang. Build PnetCDF from source (C-only):

```bash
# Download pnetcdf-1.14.1 from https://parallel-netcdf.github.io/
cd /tmp/pnetcdf-1.14.1
export OMPI_FC=flang
export OMPI_CC=clang
export OMPI_CXX=clang++

./configure --prefix=$HOME/software \
  --disable-cxx \
  --disable-shared \
  --disable-fortran \
  MPICC=mpicc \
  CFLAGS="-O3"

make -j8 && make install
```

#### PIO Library

Build PIO with LLVM compilers (see PIO documentation). Install to `$HOME/software`.

#### Building MPAS

```bash
export PKG_CONFIG_PATH="$HOME/software/lib/pkgconfig:$PKG_CONFIG_PATH"

make -j8 llvm \
  CORE=atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double
```

For `init_atmosphere`:
```bash
make -j8 llvm \
  CORE=init_atmosphere \
  PIO=$HOME/software \
  NETCDF=/opt/homebrew \
  PNETCDF=$HOME/software \
  PRECISION=double
```

#### LLVM Build Technical Details

The `llvm` target configures:
- **MPI Wrappers:** Sets `OMPI_FC=flang`, `OMPI_CC=clang`, `OMPI_CXX=clang++`
- **Fortran Flags:** Uses gfortran-compatible flags (`-fconvert=big-endian`, `-ffree-form`, `-fdefault-real-8`)
- **MPI Module:** Uses `-DNOMPIMOD` to use `include 'mpif.h'` instead of gfortran's binary `mpi.mod`

| NVIDIA/PGI Flag | LLVM flang Equivalent |
|-----------------|----------------------|
| `-Mbyteswapio` | `-fconvert=big-endian` |
| `-Mfreeform` | `-ffree-form` |
| `-r8` | `-fdefault-real-8` |
| `-Mstandard` | (not needed) |

#### Known Warnings

These warnings are expected and harmless:
- `-Wl,-flat_namespace: 'linker' input unused` - MPI wrapper passes unused linker flag
- `NUMERIC_STORAGE_SIZE from ISO_FORTRAN_ENV is not well-defined` - Due to `-fdefault-real-8`
- `Unrecognized compiler directive was ignored [-Wignored-directive]` - Intel `!DIR$` directives
- `Reference to procedure has an implicit interface` - Due to `include 'mpif.h'` usage

### Environment Variables

The Makefile reads several environment variables:

- `NETCDF` - Path to NetCDF installation
- `PNETCDF` - Path to Parallel NetCDF
- `PIO` - Path to PIO library
- `GPTL` - Path to GPTL profiling library

## Build Pipeline

The CMake build follows this sequence:

1. **Framework compilation** - `MPAS::framework`
2. **Operators compilation** - `MPAS::operators`
3. **Tool generation** - Registry parser, namelist generator
4. **External libraries** - ezxml, SMIOL, ESMF time
5. **Per-core compilation** - Atmosphere, ocean, etc.
6. **Physics libraries** - NoahMP, WRF physics, MMM physics, UGWP
7. **Executable linking** - Merge all targets
8. **Data file installation** - Physics tables, coefficient files

## Dependencies

### Required

| Dependency | Purpose |
|------------|---------|
| MPI | Distributed memory parallelism |
| PnetCDF | Parallel NetCDF I/O |

### Optional

| Dependency | Purpose |
|------------|---------|
| PIO | Parallel I/O abstraction |
| NetCDF | Additional I/O format |
| GPTL | Performance profiling |
| ESMF | Earth System Modeling Framework |

## Physics Data Files

Physics parameterizations require lookup tables and coefficient files. The build system can automatically fetch these from the MPAS-Data repository:

```bash
# Data is fetched from:
# https://github.com/MPAS-Dev/MPAS-Data (v8.0 tag)
```

Data categories:
- Microphysics lookup tables (Thompson, WSM6)
- Radiation coefficient files (RRTMG)
- Land surface parameters (Noah, NoahMP)
- Ozone climatology
- Aerosol optical properties

## CMake Modules

Custom CMake modules are in `cmake/`:

- `FindNetCDF.cmake` - Locate NetCDF
- `FindPnetCDF.cmake` - Locate Parallel NetCDF
- `FindPIO.cmake` - Locate PIO
- `MPASFunctions.cmake` - MPAS-specific CMake functions

## Troubleshooting

### Common Issues

**MPI not found:**
```bash
export MPI_HOME=/path/to/mpi
cmake -DMPI_Fortran_COMPILER=mpif90 ..
```

**NetCDF not found:**
```bash
export NETCDF=/path/to/netcdf
cmake -DNETCDF_ROOT=/path/to/netcdf ..
```

**PnetCDF not found:**
```bash
export PNETCDF=/path/to/pnetcdf
cmake -DPNETCDF_ROOT=/path/to/pnetcdf ..
```

### Verbose Build

```bash
make VERBOSE=1
```

### Clean Build

```bash
rm -rf build
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Cross-Compilation

For HPC systems with compute nodes different from login nodes:

```bash
cmake -DCMAKE_SYSTEM_NAME=Linux \
      -DCMAKE_Fortran_COMPILER=ftn \
      -DCMAKE_C_COMPILER=cc \
      ..
```

## Related Documentation

- `INSTALL` - Installation instructions
- `README.md` - Project overview
