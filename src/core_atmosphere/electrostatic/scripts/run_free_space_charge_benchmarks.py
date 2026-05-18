#!/usr/bin/env python3
"""
Free-space Gaussian charge benchmark runner for MPAS electrostatics.

This task-level runner prepares the zero-duration MPAS solve for smooth
Gaussian free-space charge sources. Error analysis and diagnostic plotting are
added by the follow-on analysis task.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import sys


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import free_space_charge_analytics as analytic  # noqa: E402
import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402
from run_tripole_supercell import ensure_output_stream_list, seed_run_dir_from_template  # noqa: E402


SOURCES = analytic.SOURCES

ELECTROSTATIC_CONFIG = (
    ("config_electrostatic_enable", ".true."),
    ("config_electrostatic_solve_at_init", ".true."),
    ("config_electrostatic_source", "'gaussian_dipole_y'"),
    ("config_poisson_preconditioner", "'jacobi'"),
    ("config_poisson_tol", "1.0e-10"),
    ("config_poisson_max_iter", "5000"),
    ("config_electrostatic_bc_ground", "0.0"),
    ("config_electrostatic_source_charge", "20.0"),
    ("config_electrostatic_source_sigma", "2000.0"),
    ("config_electrostatic_dipole_separation", "8000.0"),
)

STALE_PATTERNS = (
    "output.nc",
    "run.out",
    "log.atmosphere.*.out",
    "log.atmosphere.*.err",
)


def configure_namelist_text(
    text,
    *,
    source,
    charge,
    sigma,
    separation,
    poisson_tol,
    poisson_max_iter,
):
    """Return namelist text configured for a zero-duration Gaussian solve."""
    if source not in SOURCES:
        raise ValueError(f"unsupported Gaussian source: {source}")

    text = tier_a1.set_namelist_value(
        text,
        "nhyd_model",
        "config_run_duration",
        "'00_00:00:00'",
    )
    entries = tuple(
        (key, f"'{source}'")
        if key == "config_electrostatic_source"
        else (key, f"{poisson_tol:.1e}")
        if key == "config_poisson_tol"
        else (key, str(poisson_max_iter))
        if key == "config_poisson_max_iter"
        else (key, str(charge))
        if key == "config_electrostatic_source_charge"
        else (key, str(sigma))
        if key == "config_electrostatic_source_sigma"
        else (key, str(separation))
        if key == "config_electrostatic_dipole_separation"
        else (key, value)
        for key, value in ELECTROSTATIC_CONFIG
    )
    text = tier_a1.replace_namelist_block(text, "electrostatic", entries)
    text = re.sub(r"\n{2,}(&electrostatic)", r"\n\1", text)
    return text.rstrip() + "\n"


def configure_namelist(path, **kwargs):
    """Configure namelist.atmosphere in place."""
    path.write_text(configure_namelist_text(path.read_text(), **kwargs))


def prepare_run_dir(
    template_run_dir,
    run_dir,
    model,
    ranks,
    source,
    charge,
    sigma,
    separation,
    poisson_tol,
    poisson_max_iter,
):
    """Create and configure the isolated free-space benchmark run directory."""
    run_dir = run_dir.expanduser().resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    seed_run_dir_from_template(template_run_dir, run_dir)
    tier_a1.symlink_model(model.expanduser(), run_dir)

    namelist = run_dir / "namelist.atmosphere"
    streams = run_dir / "streams.atmosphere"
    if not namelist.exists() or not streams.exists():
        raise FileNotFoundError(
            f"{run_dir} needs namelist.atmosphere and streams.atmosphere; "
            f"seed it from a complete supercell run directory"
        )

    configure_namelist(
        namelist,
        source=source,
        charge=charge,
        sigma=sigma,
        separation=separation,
        poisson_tol=poisson_tol,
        poisson_max_iter=poisson_max_iter,
    )
    tier_a1.ensure_streams_netcdf(streams)
    ensure_output_stream_list(run_dir / "stream_list.atmosphere.output")

    partition = tier_a1.find_partition_file(run_dir, ranks, namelist.read_text())
    return run_dir, partition


def remove_stale_outputs(run_dir):
    """Remove stale output products before rerunning MPAS."""
    for pattern in STALE_PATTERNS:
        for path in run_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)


def run_mpas(run_dir, ranks, mpiexec):
    """Run MPAS in the run directory and capture stdout/stderr in run.out."""
    remove_stale_outputs(run_dir)
    command = [mpiexec, "-n", str(ranks), "./atmosphere_model"]
    with (run_dir / "run.out").open("w") as log:
        result = subprocess.run(
            command,
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"MPAS failed in {run_dir}; see {run_dir / 'run.out'}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        choices=SOURCES,
        default="gaussian_dipole_y",
        help="Gaussian charge source mode to write to config_electrostatic_source.",
    )
    parser.add_argument(
        "--template-run-dir",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/supercell"),
        help="Existing supercell run directory to seed inputs from.",
    )
    parser.add_argument(
        "--run-dir",
        type=pathlib.Path,
        default=None,
        help="Isolated run directory to create or reuse.",
    )
    parser.add_argument(
        "--model",
        type=pathlib.Path,
        default=pathlib.Path("./atmosphere_model"),
        help="Path to the atmosphere_model executable.",
    )
    parser.add_argument("--ranks", type=int, default=8, help="MPI ranks to use.")
    parser.add_argument("--mpiexec", default="mpiexec", help="MPI launcher executable.")
    parser.add_argument(
        "--charge",
        type=float,
        default=20.0,
        help="Gaussian lobe charge written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=2000.0,
        help="Gaussian width in meters written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--separation",
        type=float,
        default=8000.0,
        help="Dipole lobe separation in meters written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--boundary-margin",
        type=float,
        default=12000.0,
        help="Interior comparison margin in meters for the follow-on analysis step.",
    )
    parser.add_argument(
        "--core-radius",
        type=float,
        default=None,
        help="Optional source-core exclusion radius in meters for follow-on analysis.",
    )
    parser.add_argument(
        "--e-relative-floor",
        type=float,
        default=0.0,
        help="Relative-error floor for electric-field analysis.",
    )
    parser.add_argument(
        "--poisson-tol",
        type=float,
        default=1.0e-10,
        help="PCG tolerance written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--poisson-max-iter",
        type=int,
        default=5000,
        help="PCG iteration cap written to namelist.atmosphere.",
    )
    parser.add_argument(
        "--summary",
        type=pathlib.Path,
        default=None,
        help="Output JSON summary path for the follow-on analysis step.",
    )
    parser.add_argument(
        "--plot",
        type=pathlib.Path,
        default=None,
        help="Output PNG diagnostics path for the follow-on analysis step.",
    )
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Skip MPAS execution and analyze an existing output.nc.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare the run directory but do not run or analyze output.",
    )
    args = parser.parse_args(argv)

    root = pathlib.Path("~/Data/MPAS/poisson_free_space_charge")
    if args.run_dir is None:
        args.run_dir = root / args.source / "run"
    if args.summary is None:
        args.summary = root / "results" / f"{args.source}_summary.json"
    if args.plot is None:
        args.plot = root / "results" / f"{args.source}_diagnostics.png"
    return args
