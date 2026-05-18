#!/usr/bin/env python3
"""
Prepare and run a diagnostic charge-coupled supercell case.

This helper enables the Phase 1 one-way charge coupling path:
MPAS state -> stub rho_charge -> Poisson solve -> phi/E diagnostics.
It does not feed electric fields back into dynamics or microphysics.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402
from run_tripole_supercell import seed_run_dir_from_template  # noqa: E402


ELECTROSTATIC_CONFIG = (
    ("config_electrostatic_enable", ".true."),
    ("config_electrostatic_solve_at_init", ".false."),
    ("config_electrostatic_source", "'stub'"),
    ("config_electrostatic_interval", "60.0"),
    ("config_stub_alpha", "1.0e-8"),
    ("config_stub_beta", "1.0e-8"),
    ("config_poisson_preconditioner", "'jacobi'"),
    ("config_poisson_tol", "1.0e-10"),
    ("config_poisson_max_iter", "5000"),
    ("config_electrostatic_bc_ground", "0.0"),
)

COUPLED_OUTPUT_FIELDS = (
    "xtime",
    "xCell",
    "yCell",
    "zgrid",
    "areaCell",
    "w",
    "scalars",
    "rho_charge",
    "phi",
    "E_normal",
    "E_vector",
    "cg_iter_count",
    "cg_residual_initial",
    "cg_residual_final",
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
    run_duration,
    interval,
    stub_alpha,
    stub_beta,
    poisson_tol,
    poisson_max_iter,
    solve_at_init,
):
    """Return namelist text configured for the dynamic stub source."""
    text = tier_a1.set_namelist_value(
        text,
        "nhyd_model",
        "config_run_duration",
        f"'{run_duration}'",
    )
    entries = tuple(
        (key, ".true." if solve_at_init else ".false.")
        if key == "config_electrostatic_solve_at_init"
        else (key, str(interval))
        if key == "config_electrostatic_interval"
        else (key, str(stub_alpha))
        if key == "config_stub_alpha"
        else (key, str(stub_beta))
        if key == "config_stub_beta"
        else (key, str(poisson_tol))
        if key == "config_poisson_tol"
        else (key, str(poisson_max_iter))
        if key == "config_poisson_max_iter"
        else (key, value)
        for key, value in ELECTROSTATIC_CONFIG
    )
    text = tier_a1.replace_namelist_block(text, "electrostatic", entries)
    text = re.sub(r"\n{2,}(&electrostatic)", r"\n\1", text)
    return text.rstrip() + "\n"


def configure_namelist(path, **kwargs):
    """Configure namelist.atmosphere in place."""
    path.write_text(configure_namelist_text(path.read_text(), **kwargs))


def ensure_coupled_output_stream_list(stream_list_path):
    """Write a compact output stream list for dynamic coupling inspection."""
    present = set()
    lines = []
    for field in COUPLED_OUTPUT_FIELDS:
        if field not in present:
            lines.append(field)
            present.add(field)

    stream_list_path.write_text("\n".join(lines).rstrip() + "\n")


def configure_streams(streams_path, output_interval):
    """Force netCDF streams and write dynamic samples to output.nc."""
    tier_a1.ensure_streams_netcdf(streams_path)
    tree = ET.parse(streams_path)
    root = tree.getroot()
    output = root.find("./stream[@name='output']")
    if output is None:
        output = ET.SubElement(root, "stream", {"name": "output", "type": "output"})

    output.set("type", "output")
    output.set("io_type", "netcdf")
    output.set("filename_template", "output.nc")
    output.set("filename_interval", "none")
    output.set("output_interval", output_interval)

    if output.find("./file[@name='stream_list.atmosphere.output']") is None:
        ET.SubElement(output, "file", {"name": "stream_list.atmosphere.output"})

    ET.indent(tree, space="  ")
    tree.write(streams_path, encoding="unicode")


def prepare_run_dir(
    *,
    template_run_dir,
    run_dir,
    model,
    ranks,
    run_duration,
    interval,
    stub_alpha,
    stub_beta,
    poisson_tol,
    poisson_max_iter,
    solve_at_init,
    output_interval,
):
    """Create and configure an isolated dynamic charge-coupled run directory."""
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
        run_duration=run_duration,
        interval=interval,
        stub_alpha=stub_alpha,
        stub_beta=stub_beta,
        poisson_tol=poisson_tol,
        poisson_max_iter=poisson_max_iter,
        solve_at_init=solve_at_init,
    )
    configure_streams(streams, output_interval)
    ensure_coupled_output_stream_list(run_dir / "stream_list.atmosphere.output")
    partition = tier_a1.find_partition_file(run_dir, ranks, namelist.read_text())
    return run_dir, partition


def remove_stale_outputs(run_dir):
    """Remove stale MPAS products before rerunning."""
    for pattern in STALE_PATTERNS:
        for path in run_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)


def run_mpas(run_dir, ranks, mpiexec):
    """Run MPAS and capture combined stdout/stderr in run.out."""
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
        "--template-run-dir",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/supercell"),
        help="Existing supercell run directory to seed inputs from.",
    )
    parser.add_argument(
        "--run-dir",
        type=pathlib.Path,
        default=pathlib.Path("~/Data/MPAS/poisson_charge_coupled_supercell/run"),
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
        "--run-duration",
        default="00_00:10:00",
        help="MPAS config_run_duration value.",
    )
    parser.add_argument(
        "--electrostatic-interval",
        type=float,
        default=60.0,
        help="Seconds between diagnostic electrostatic solves.",
    )
    parser.add_argument("--stub-alpha", type=float, default=1.0e-8)
    parser.add_argument("--stub-beta", type=float, default=1.0e-8)
    parser.add_argument("--poisson-tol", type=float, default=1.0e-10)
    parser.add_argument("--poisson-max-iter", type=int, default=5000)
    parser.add_argument(
        "--output-interval",
        default="00:10:00",
        help="MPAS output stream interval for output.nc.",
    )
    parser.add_argument(
        "--solve-at-init",
        action="store_true",
        help="Also run a diagnostic electrostatic solve during initialization.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare the run directory but do not run MPAS.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run_dir, partition = prepare_run_dir(
        template_run_dir=args.template_run_dir,
        run_dir=args.run_dir,
        model=args.model,
        ranks=args.ranks,
        run_duration=args.run_duration,
        interval=args.electrostatic_interval,
        stub_alpha=args.stub_alpha,
        stub_beta=args.stub_beta,
        poisson_tol=args.poisson_tol,
        poisson_max_iter=args.poisson_max_iter,
        solve_at_init=args.solve_at_init,
        output_interval=args.output_interval,
    )
    print(f"Prepared {run_dir}")
    print(f"Using partition {partition}")
    if args.prepare_only:
        return 0
    run_mpas(run_dir, args.ranks, args.mpiexec)
    print(f"Wrote {run_dir / 'output.nc'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
