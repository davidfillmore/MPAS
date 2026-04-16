# MPAS

Clean fork of [MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model).
**No chemistry** — this is the baseline MPAS, distinct from the sibling
`DAVINCI-MPAS` fork (DAVINCI chemistry, in-tree Fortran) and
`MPAS-Model-ACOM-dev` fork (MUSICA/MICM chemistry integration).

## ON SESSION START

**The first thing to do in any new session in this repo is read
[`../MPAS-Papers/CLAUDE.md`](../MPAS-Papers/CLAUDE.md).** That file holds the
shared workflow rules, planning conventions, ecosystem overview, and the
centralized history/journal practice for all MPAS-ecosystem work. Do this
before answering any question or touching code.

## Project Guidance

Plans, papers, and session history for this repo live in
[`../MPAS-Papers/`](../MPAS-Papers/) — the private planning hub. Nothing about
planning or journaling is kept in this public code fork.

Key workflow rules (full detail in `../MPAS-Papers/CLAUDE.md`):

- **Never begin implementation without explicit instructions to proceed** —
  present plans and wait.
- **Stop for discussion between plan phases** — report results after each
  phase (build, run, validate) before starting the next.
- **Journal at end of session** into `../MPAS-Papers/journal/YYYY-MM-DD.md`
  and update `../MPAS-Papers/HISTORY.md`.

## Project History

Session history is **not** kept in this repo. Journal entries for work done
here land in [`../MPAS-Papers/journal/YYYY-MM-DD.md`](../MPAS-Papers/journal/)
and are indexed from [`../MPAS-Papers/HISTORY.md`](../MPAS-Papers/HISTORY.md).
See [`../MPAS-Papers/CLAUDE.md`](../MPAS-Papers/CLAUDE.md) § "Project History"
for the full convention. Session logs (for the AI Compute Log) for sessions run
from this repo are under `~/.claude/projects/-Users-fillmore-EarthSystem-MPAS/`.

## Build and Run

| Doc | Purpose |
|-----|---------|
| [BUILD.md](BUILD.md) | CMake and legacy Makefile builds (GNU/Ubuntu conda-forge, LLVM/macOS Homebrew) |
| [RUN.md](RUN.md) | Test cases: supercell, mountain wave, Jablonowski-Williamson baroclinic wave |

## Repository Setup

- **origin**: `ssh://git@github.com/davidfillmore/MPAS` (fork)
- **upstream**: `https://github.com/MPAS-Dev/MPAS-Model` (parent)

### Branches

- `master` — tracks upstream master (stable releases)
- `develop` — tracks upstream develop (active development; current working branch)

### Syncing with Upstream

```sh
git checkout develop
git fetch upstream develop
git merge upstream/develop
git push origin develop
```

Same pattern for `master`.

## Host-Specific Notes

**macOS (Homebrew LLVM).** PIO on this host was built without PnetCDF, so every
stream in `streams.atmosphere` and `streams.init_atmosphere` must set
`io_type="netcdf"` explicitly for multi-rank runs — otherwise PIO picks a
parallel default, finds no usable backend, and aborts at stream open. See
[RUN.md](RUN.md) § "Important: I/O Configuration".

**Ubuntu (conda-forge gfortran).** All compilers and libraries (except PIO)
come from the `mpas` conda environment. Activate before building:
`conda activate mpas`.

## Key Source Locations

| Component | Path |
|-----------|------|
| Framework | `src/framework/` |
| Operators | `src/operators/` |
| Atmosphere core | `src/core_atmosphere/` |
| Init atmosphere | `src/core_init_atmosphere/` |
| Dynamics | `src/core_atmosphere/dynamics/` |
| Physics | `src/core_atmosphere/physics/` |
| Registry (metadata) | `src/core_atmosphere/Registry.xml` |

## Run Directory

Test cases are run from `~/Data/MPAS/`, separate from the source tree:

| Directory | Description |
|-----------|-------------|
| `~/Data/MPAS/supercell` | Idealized supercell thunderstorm (Kessler microphysics) |
| `~/Data/MPAS/mountain_wave` | 2D Schaer mountain wave |
| `~/Data/MPAS/jw_baroclinic_wave` | Global Jablonowski-Williamson baroclinic wave |

## Notes for AI Assistants

1. **Fortran 2008 codebase, MPI parallelism.** Registry-driven (`Registry.xml`
   parsed at build time).
2. **Never mix compilers.** flang `.mod` files are incompatible with gfortran
   `.mod` files.
3. **Symlink the executable** from the run directory to the repo's
   `atmosphere_model` — avoid `cp` so rebuilds are picked up automatically.
4. **Clobber mode matters.** MPAS defaults to `clobber_mode = never_modify`. If
   `output.nc` already exists, the model silently skips all output writes. Move
   or remove it before re-running.
5. **LANDUSE.TBL is required** in the run directory even when
   `config_physics_suite='none'`.
