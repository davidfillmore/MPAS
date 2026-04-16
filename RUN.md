# Running MPAS Test Cases

This document describes how to run MPAS atmosphere test cases.

## Supercell Test Case

The supercell thunderstorm is an idealized convection test case located at `~/Data/MPAS/supercell`.

### Important: I/O Configuration

The `streams.atmosphere` file must use `io_type="netcdf"` for input and output streams to avoid PnetCDF compatibility issues:

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

### Prerequisites

1. Build the atmosphere model (see [BUILD.md](BUILD.md)):
   ```bash
   make -j8 gfortran \
     CORE=atmosphere \
     PIO="$PIO" NETCDF="$NETCDF" PNETCDF="$PNETCDF" PRECISION=double
   ```

2. Verify the executable exists:
   ```bash
   ls -la atmosphere_model
   ```

### Running the Test

**Important:** You must remove or move any existing `output.nc` before running. MPAS defaults to `clobber_mode = never_modify`, so if `output.nc` already exists the model will silently skip all output writes — the run completes but produces no new data.

```bash
cd ~/Data/MPAS/supercell

# Archive previous run output (REQUIRED — model won't overwrite existing output.nc)
timestamp=$(date +%Y%m%d_%H%M%S)
[ -f output.nc ] && mv output.nc output.${timestamp}.nc
[ -f log.atmosphere.0000.out ] && mv log.atmosphere.0000.out log.atmosphere.0000.${timestamp}.out

# Run with 8 MPI ranks (recommended for 10-core machine)
mpiexec -n 8 ~/EarthSystem/MPAS/atmosphere_model
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

### Observed behavior on this host

The following values were captured from a clean 8-rank run of the
supercell case on this machine (conda-forge gfortran build, see
[BUILD.md](BUILD.md)). Use them as a sanity check that a fresh build is
producing the right answer.

| Metric | Value |
|--------|-------|
| Ranks | 8 (`supercell.graph.info.part.8`) |
| Timestep | 3.0 s |
| Total steps | 2400 (for the default 2-hour run) |
| Wall time per integration step | ≈ 0.47 s |
| Full-run wall time (projected) | ≈ 19 min |
| `output.nc` size at completion | ≈ 4.6 GB |

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

## Test Case Data

Idealized test case data is downloaded from
[NCAR MPAS v7.0](https://www2.mmm.ucar.edu/projects/mpas/site/access_code/idealized.html).
