# Running MPAS Test Cases

This document describes how to run MPAS atmosphere test cases.

## Supercell Test Case

The supercell thunderstorm is an idealized convection test case located at `~/Data/MPAS/supercell`.

### Important: I/O Configuration

On hosts where PIO was built without PnetCDF (e.g. the macOS LLVM build
described in [BUILD.md](BUILD.md)), every input and output stream in
**both** `streams.init_atmosphere` and `streams.atmosphere` must set
`io_type="netcdf"` explicitly. Otherwise PIO picks a parallel default,
finds no usable backend, and aborts with
`CRITICAL ERROR: Could not open input file ...`.

```xml
<immutable_stream name="input"
                  type="input"
                  io_type="netcdf"
                  filename_template="supercell_init.nc"
                  input_interval="initial_only"/>

<stream name="output"
        type="output"
        io_type="netcdf"
        ...>
```

On hosts where PIO does have PnetCDF support, the explicit `io_type` is
harmless — it just forces the serial path. Leave it in for portability.

### Prerequisites

1. Build both the init_atmosphere and atmosphere cores (see
   [BUILD.md](BUILD.md)). The legacy Makefile builds one core per
   invocation:
   ```bash
   make -j8 gfortran CORE=init_atmosphere \
     PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
   make -j8 gfortran CORE=atmosphere \
     PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
   ```

   On macOS, substitute `llvm` for `gfortran` and drop `PNETCDF=` if PIO
   was built without PnetCDF.

2. Verify the executables exist:
   ```bash
   ls -la init_atmosphere_model atmosphere_model
   ```

3. `LANDUSE.TBL` must be present in the run directory. Physics init
   reads it even when `config_physics_suite='none'`; without it the run
   aborts with `subroutine landuse_init_forMPAS: failure opening LANDUSE.TBL`.

### Running the Test

The supercell case runs in two stages: `init_atmosphere_model` generates
initial conditions from the mesh, then `atmosphere_model` integrates
forward.

```bash
cd ~/Data/MPAS/supercell

# Stage 1: initial conditions (produces supercell_init.nc)
rm -f log.init_atmosphere.*.{out,err} supercell_init.nc
mpiexec -n 8 ./init_atmosphere_model
```

**Important:** You must remove or move any existing `output.nc` before
running stage 2. MPAS defaults to `clobber_mode = never_modify`, so if
`output.nc` already exists the model will silently skip all output
writes — the run completes but produces no new data.

```bash
# Stage 2: atmosphere forecast
timestamp=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${timestamp}.nc
[ -f log.atmosphere.0000.out ] && mv log.atmosphere.0000.out log.atmosphere.0000.${timestamp}.out

# Run with 8 MPI ranks (recommended for 10-core machine)
mpiexec -n 8 ./atmosphere_model
```

### MPI Rank Selection

Choose based on available cores and partition files:

| Ranks | Partition File | Notes |
|-------|----------------|-------|
| 8 | `supercell.graph.info.part.8` | Recommended for 10-core machine |
| 12 | `supercell.graph.info.part.12` | Slight oversubscription |
| 16 | `supercell.graph.info.part.16` | Oversubscribed on 10 cores |

The partition file prefix is set in `namelist.atmosphere`:
```
config_block_decomp_file_prefix = 'supercell.graph.info.part.'
```

### Configuration

Key settings in `namelist.atmosphere`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `config_dt` | 3.0 | Timestep (seconds) |
| `config_run_duration` | '02:00:00' | Run length (2 hours) |
| `config_microp_scheme` | 'mp_kessler' | Microphysics scheme |

### Output Files

| File | Description |
|------|-------------|
| `output.nc` | Model output (written every 2 minutes) |
| `log.atmosphere.0000.out` | Run log with diagnostics |

### Verifying the Run

Check the log for successful completion:
```bash
tail -20 log.atmosphere.0000.out
```

Look for final timestep and timing statistics.

### Observed behavior

Reference numbers from clean 8-rank runs of the supercell case. Use them
as a sanity check that a fresh build is producing the right answer.

**Ubuntu host, conda-forge gfortran build** (2-hour run):

| Metric | Value |
|--------|-------|
| Ranks | 8 (`supercell.graph.info.part.8`) |
| Timestep | 3.0 s |
| Total steps | 2400 |
| Wall time per integration step | ≈ 0.47 s |
| Full-run wall time (projected) | ≈ 19 min |
| `output.nc` size at completion | ≈ 4.6 GB |

**macOS host, Homebrew LLVM flang build** (3-minute smoke-test run):

| Metric | Value |
|--------|-------|
| Ranks | 8 (`supercell.graph.info.part.8`) |
| Timestep | 3.0 s |
| Total steps | 60 |
| `init_atmosphere` wall time | ≈ 2 s |
| `atmosphere` wall time | ≈ 23 s |
| `supercell_init.nc` size | ≈ 357 MB |
| `output.nc` size (2 Time records) | ≈ 517 MB |

During the first ~50 minutes of simulated time, typical per-step
diagnostics look like:

```
global min, max w       -20.75 ...  47.60   (m/s)
global min, max u       -41.42 ...  41.43   (m/s)
global min, max scalar 1  0.0  ...  0.0157  (qv)
global min, max scalar 2  0.0  ...  0.0049  (qc)
global min, max scalar 3  0.0  ...  0.0176  (qr)
Timing for integration step: ~0.47 s
```

A strong updraft (`w` peaking around +47 m/s by t ≈ 50 min) and the
Kessler scalars (`qv`, `qc`, `qr`) all non-negative are the two quick
signs that the dynamics-transport coupling and microphysics are working.

The typical smoke test used in this repo is to run through 50–60
simulated minutes (≈ 7–8 min wall time) and confirm these ranges, then
kill the run. A full 2-hour simulation is only needed when a complete
`output.nc` is required.

### Quick Run Script

Create a helper script `run.sh`:
```bash
#!/bin/bash
cd ~/Data/MPAS/supercell

# Archive previous output
ts=$(date +%Y%m%d_%H%M%S)
for f in output.nc log.atmosphere.*.out; do
    [ -f "$f" ] && mv "$f" "${f%.nc}.${ts}.nc" 2>/dev/null || mv "$f" "${f%.out}.${ts}.out" 2>/dev/null
done

# Run
mpiexec -n 8 ~/EarthSystem/MPAS/atmosphere_model 2>&1 | tee run.out
```

## Mountain Wave Test Case

A 2D mountain wave (Schaer test) in `~/Data/MPAS/mountain_wave/`. Tests
non-hydrostatic dynamics over orography with 70 vertical levels.

### Configuration

| Parameter | Value |
|-----------|-------|
| `config_dt` | 6.0 s |
| `config_run_duration` | 5 hours |
| `config_nvertlevels` | 70 |
| `config_physics_suite` | none |
| `config_scalar_advection` | false |

### Running

```bash
cd ~/Data/MPAS/mountain_wave
rm -f output.nc log.atmosphere.*.out
mpiexec -n 4 ./atmosphere_model 2>&1 | tee run.out
```

Partition files available for 2, 4, 6, 8 ranks.

## Jablonowski-Williamson Baroclinic Wave

A global baroclinic instability test on a 120-km (40,962 cell) mesh in
`~/Data/MPAS/jw_baroclinic_wave/`. This is a standard dynamical core
benchmark (Jablonowski & Williamson 2006).

### Configuration

| Parameter | Value |
|-----------|-------|
| `config_dt` | 450.0 s |
| `config_run_duration` | 16 days |
| `config_nvertlevels` | 26 |
| `config_horiz_mixing` | 2d_smagorinsky |
| `config_physics_suite` | none |
| `config_scalar_advection` | false |

**Note:** The default 16-day run is computationally expensive. For a quick
verification, reduce `config_run_duration` to `'4_00:00:00'` (4 days).

### Running

```bash
cd ~/Data/MPAS/jw_baroclinic_wave
rm -f output.nc log.atmosphere.*.out
mpiexec -n 8 ./atmosphere_model 2>&1 | tee run.out
```

Partition files available for 2, 4, 6, 8, 12, 16, 24 ranks.

## Poisson Tier A.1 MMS

The electrostatic Phase 1A benchmark runner prepares Cartesian MMS run
directories under `~/Data/MPAS/poisson_tier_A1`, runs `atmosphere_model`,
and writes `results/tier_A1_convergence.csv` plus
`results/tier_A1_convergence.png`.

```bash
~/miniconda3/envs/mpas/bin/python \
  src/core_atmosphere/electrostatic/scripts/run_tier_A1_cartesian_mms.py \
  --run-root ~/Data/MPAS/poisson_tier_A1 \
  --mesh-list 15km 7.5km 3.75km 1.875km \
  --ranks 4 \
  --model ./atmosphere_model
```

Each mesh needs either `meshes/<mesh>/` or `runs/<mesh>/` to contain
`namelist.atmosphere`, `streams.atmosphere`, the input NetCDF files, and
a partition file matching `--ranks`, for example
`*.graph.info.part.4` when using `--ranks 4`. The runner symlinks the
model executable, forces stream `io_type="netcdf"`, removes stale
`output.nc` and MPAS logs before reruns, and writes all generated run
artifacts outside the source tree.

## Test Case Data

Idealized test case data is downloaded from
[NCAR MPAS v7.0](https://www2.mmm.ucar.edu/projects/mpas/site/access_code/idealized.html).
