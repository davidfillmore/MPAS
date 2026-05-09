#!/usr/bin/env python3
"""
Generate Tier A.2 spherical MMS mesh/input bundles.

Each bundle under run-root/meshes/<mesh> contains:

* grid.nc
* graph.info and graph.info.part.<ranks>
* init.nc from init_atmosphere_model
* namelist/streams files used by run_tier_A2_sphere_mms.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys

import netCDF4 as nc
import numpy as np


EARTH_RADIUS_M = 6371229.0
DEFAULT_MESHES = ("960km", "480km", "240km")
OUTPUT_FIELDS = (
    "latCell",
    "lonCell",
    "zgrid",
    "areaCell",
    "phi",
    "rho_charge",
    "E_normal",
    "E_vector",
    "cg_iter_count",
    "cg_residual_initial",
    "cg_residual_final",
)


def mesh_spacing_km(mesh_name):
    """Infer mesh spacing in kilometers from labels such as 240km."""
    name = mesh_name.lower().replace("_", ".")
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*km", name)
    if match:
        return float(match.group(1))
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*m", name)
    if match:
        return float(match.group(1)) / 1000.0
    raise ValueError(f"Cannot infer mesh spacing from {mesh_name!r}")


def write_text_if_needed(path, text, force):
    """Write a text file unless it already exists and force is false."""
    if force or not path.exists():
        path.write_text(text.rstrip() + "\n")


def init_namelist_text(nvertlevels, ztop):
    """Return a minimal init_atmosphere namelist for vertical-grid generation."""
    return f"""
&nhyd_model
    config_init_case = 5
    config_start_time = '0000-01-01_00:00:00'
    config_stop_time = '0000-01-01_00:00:00'
    config_theta_adv_order = 3
    config_coef_3rd_order = 0.25
    config_interface_projection = 'linear_interpolation'
/
&dimensions
    config_nvertlevels = {nvertlevels}
    config_nsoillevels = 4
    config_nfglevels = 38
    config_nfgsoillevels = 4
    config_gocartlevels = 30
/
&data_sources
    config_geog_data_path = './'
/
&vertical_grid
    config_ztop = {ztop:.1f}
    config_hybrid_coordinate = false
    config_hybrid_top_z = {ztop:.1f}
    config_nsmterrain = 1
    config_smooth_surfaces = false
    config_dzmin = 0.3
    config_nsm = 30
    config_tc_vertical_grid = false
    config_blend_bdy_terrain = false
/
&preproc_stages
    config_static_interp = false
    config_native_gwd_static = false
    config_native_gwd_gsl_static = false
    config_vertical_grid = true
    config_met_interp = false
    config_input_sst = false
    config_frac_seaice = false
/
&io
    config_pio_num_iotasks = 0
    config_pio_stride = 1
/
&decomposition
    config_block_decomp_file_prefix = 'graph.info.part.'
/
"""


def init_streams_text():
    """Return streams.init_atmosphere for the generated mesh bundle."""
    return """
<streams>
  <immutable_stream name="input" type="input" io_type="netcdf" filename_template="grid.nc" input_interval="initial_only" />
  <immutable_stream name="output" type="output" io_type="netcdf" filename_template="init.nc" packages="initial_conds" output_interval="initial_only" />
</streams>
"""


def atmosphere_namelist_text():
    """Return a zero-duration MPAS-A namelist for the electrostatic MMS solve."""
    return """
&nhyd_model
    config_dt = 3.0
    config_start_time = '0000-01-01_00:00:00'
    config_split_dynamics_transport = false
    config_number_of_sub_steps = 6
    config_dynamics_split_steps = 1
    config_h_mom_eddy_visc2  = 500.0
    config_h_mom_eddy_visc4  = 0.0
    config_v_mom_eddy_visc2  = 500.0
    config_h_theta_eddy_visc2 = 500.0
    config_h_theta_eddy_visc4 = 0.0
    config_v_theta_eddy_visc2 = 500.0
    config_horiz_mixing = '2d_fixed'
    config_len_disp = 538.86
    config_theta_adv_order  = 3
    config_w_adv_order      = 3
    config_u_vadv_order     = 3
    config_w_vadv_order     = 3
    config_theta_vadv_order = 3
    config_coef_3rd_order = 0.25
    config_epssm = 0.1
    config_smdiv = 0.1
    config_mix_full = false
    config_monotonic = true
    config_h_ScaleWithMesh = false
    config_run_duration = '00_00:00:00'
/
&damping
    config_zd = 20000.0
    config_xnutr = 0.0
/
&limited_area
    config_apply_lbcs = false
/
&io
    config_pio_num_iotasks = 0
    config_pio_stride = 1
/
&decomposition
    config_block_decomp_file_prefix = 'graph.info.part.'
/
&restart
    config_do_restart = false
/
&printout
    config_print_global_minmax_vel = true
    config_print_global_minmax_sca = true
/
&physics
    config_o3climatology = false
    config_sst_update = false
    config_sstdiurn_update = false
    config_deepsoiltemp_update = false
    config_radtlw_interval = '00:30:00'
    config_radtsw_interval = '00:30:00'
    config_bucket_update = 'none'
    config_physics_suite = 'none'
    config_microp_scheme = 'off'
/
&electrostatic
    config_electrostatic_enable = .true.
    config_electrostatic_solve_at_init = .true.
    config_electrostatic_source = 'mms_sphere_horizontal'
    config_poisson_preconditioner = 'jacobi'
    config_poisson_tol = 1.0e-8
    config_poisson_max_iter = 2000
    config_electrostatic_bc_ground = 0.0
/
"""


def atmosphere_streams_text():
    """Return streams.atmosphere for the zero-duration MMS run."""
    return """
<streams>
  <immutable_stream name="input" type="input" io_type="netcdf" filename_template="init.nc" input_interval="initial_only" />
  <immutable_stream name="restart" type="input;output" io_type="netcdf" filename_template="restart.$Y-$M-$D_$h.$m.$s.nc" input_interval="initial_only" output_interval="0:15:00" />
  <stream name="output" type="output" io_type="netcdf" filename_template="output.nc" filename_interval="none" output_interval="initial_only">
    <file name="stream_list.atmosphere.output" />
  </stream>
  <stream name="diagnostics" type="output" io_type="netcdf" filename_template="diag.$Y-$M-$D_$h.$m.$s.nc" output_interval="none" />
  <stream name="surface" type="input" io_type="netcdf" filename_template="sfc_update.nc" filename_interval="none" input_interval="none" />
  <immutable_stream name="iau" type="input" io_type="netcdf" filename_template="AmB.$Y-$M-$D_$h.$m.$s.nc" filename_interval="none" packages="iau" input_interval="initial_only" />
  <immutable_stream name="lbc_in" type="input" io_type="netcdf" filename_template="lbc.$Y-$M-$D_$h.$m.$s.nc" filename_interval="input_interval" packages="limited_area" input_interval="none" />
</streams>
"""


def output_stream_list_text():
    """Return stream_list.atmosphere.output for Tier A.2 diagnostics."""
    return "\n".join(OUTPUT_FIELDS)


def clean_generated_mesh_files(work_dir):
    """Remove stale generated mesh files from a previous failed attempt."""
    for path in work_dir.glob("*"):
        if path.is_file() or path.is_symlink():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)


def generate_spherical_mesh(bundle_dir, work_dir, spacing_km, earth_radius, plot_cell_width, force):
    """Generate grid.nc and graph.info for one uniform spherical mesh."""
    from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh

    grid = bundle_dir / "grid.nc"
    graph = bundle_dir / "graph.info"
    if grid.exists() and graph.exists() and not force:
        return

    work_dir.mkdir(parents=True, exist_ok=True)
    clean_generated_mesh_files(work_dir)

    lat = np.linspace(-90.0, 90.0, 19)
    lon = np.linspace(-180.0, 180.0, 37)
    cell_width = spacing_km * np.ones((lat.size, lon.size))

    previous = pathlib.Path.cwd()
    try:
        os.chdir(work_dir)
        build_spherical_mesh(
            cell_width,
            lon,
            lat,
            earth_radius=earth_radius,
            out_filename="base_mesh.nc",
            plot_cellWidth=plot_cell_width,
        )
    finally:
        os.chdir(previous)

    shutil.copy2(work_dir / "base_mesh.nc", grid)
    shutil.copy2(work_dir / "graph.info", graph)


def run_command(command, cwd, log_path=None):
    """Run a command, optionally teeing stdout/stderr into a log file."""
    if log_path is None:
        subprocess.run(command, cwd=cwd, check=True)
        return

    with log_path.open("w") as log:
        result = subprocess.run(
            command,
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed in {cwd}; see {log_path}")


def partition_mesh(bundle_dir, ranks, force):
    """Create a METIS partition file for one mesh."""
    part_file = bundle_dir / f"graph.info.part.{ranks}"
    if part_file.exists() and not force:
        return
    run_command(["gpmetis", "graph.info", str(ranks)], cwd=bundle_dir)


def write_bundle_inputs(bundle_dir, nvertlevels, ztop, force):
    """Write namelist and stream files used by init and atmosphere runs."""
    write_text_if_needed(
        bundle_dir / "namelist.init_atmosphere",
        init_namelist_text(nvertlevels, ztop),
        force,
    )
    write_text_if_needed(bundle_dir / "streams.init_atmosphere", init_streams_text(), force)
    write_text_if_needed(
        bundle_dir / "namelist.atmosphere", atmosphere_namelist_text(), force
    )
    write_text_if_needed(bundle_dir / "streams.atmosphere", atmosphere_streams_text(), force)
    write_text_if_needed(
        bundle_dir / "stream_list.atmosphere.output",
        output_stream_list_text(),
        force,
    )


def run_init_atmosphere(bundle_dir, init_model, ranks, mpiexec, force):
    """Run init_atmosphere_model to create init.nc."""
    output = bundle_dir / "init.nc"
    if output.exists() and not force:
        return
    for pattern in ("init.nc", "init.out", "log.init_atmosphere.*.out", "log.init_atmosphere.*.err"):
        for path in bundle_dir.glob(pattern):
            path.unlink()
    run_command(
        [mpiexec, "-n", str(ranks), str(init_model)],
        cwd=bundle_dir,
        log_path=bundle_dir / "init.out",
    )
    if not output.exists():
        raise FileNotFoundError(f"{output} was not created")


def mesh_summary(bundle_dir, mesh_name, spacing_km, ranks, nvertlevels):
    """Return summary metadata for one generated bundle."""
    with nc.Dataset(bundle_dir / "grid.nc") as dataset:
        n_cells = len(dataset.dimensions["nCells"])
        n_edges = len(dataset.dimensions["nEdges"])
        n_vertices = len(dataset.dimensions["nVertices"])
    return {
        "mesh": mesh_name,
        "spacing_km": spacing_km,
        "ranks": ranks,
        "nCells": n_cells,
        "nEdges": n_edges,
        "nVertices": n_vertices,
        "nVertLevels": nvertlevels,
        "bundle": str(bundle_dir),
    }


def setup_mesh(args, mesh_name):
    """Set up one Tier A.2 mesh bundle."""
    spacing_km = mesh_spacing_km(mesh_name)
    run_root = args.run_root.expanduser().resolve()
    bundle_dir = run_root / "meshes" / mesh_name
    work_dir = run_root / "mesh_scratch" / mesh_name
    bundle_dir.mkdir(parents=True, exist_ok=True)

    print(f"{mesh_name}: generating {spacing_km:g} km spherical mesh")
    generate_spherical_mesh(
        bundle_dir,
        work_dir,
        spacing_km,
        args.earth_radius,
        args.plot_cell_width,
        args.force,
    )
    write_bundle_inputs(bundle_dir, args.nvertlevels, args.ztop, args.force)
    print(f"{mesh_name}: partitioning for {args.ranks} ranks")
    partition_mesh(bundle_dir, args.ranks, args.force)

    if not args.mesh_only:
        init_model = args.init_model.expanduser().resolve()
        if not init_model.exists():
            raise FileNotFoundError(f"init model does not exist: {init_model}")
        print(f"{mesh_name}: running init_atmosphere_model")
        run_init_atmosphere(bundle_dir, init_model, args.ranks, args.mpiexec, args.force)

    return mesh_summary(bundle_dir, mesh_name, spacing_km, args.ranks, args.nvertlevels)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=pathlib.Path,
        required=True,
        help="Root where meshes/, runs/, results/, and mesh_scratch/ live.",
    )
    parser.add_argument(
        "--mesh-list",
        nargs="+",
        default=list(DEFAULT_MESHES),
        help="Mesh labels such as 960km 480km 240km.",
    )
    parser.add_argument("--ranks", type=int, default=8, help="MPI ranks to partition/run.")
    parser.add_argument(
        "--init-model",
        type=pathlib.Path,
        default=pathlib.Path("./init_atmosphere_model"),
        help="Path to init_atmosphere_model.",
    )
    parser.add_argument("--mpiexec", default="mpiexec", help="MPI launcher executable.")
    parser.add_argument("--nvertlevels", type=int, default=16, help="Vertical levels.")
    parser.add_argument("--ztop", type=float, default=30000.0, help="Model top in meters.")
    parser.add_argument(
        "--earth-radius",
        type=float,
        default=EARTH_RADIUS_M,
        help="Sphere radius in meters.",
    )
    parser.add_argument(
        "--mesh-only",
        action="store_true",
        help="Generate grid and partition files, but skip init.nc creation.",
    )
    parser.add_argument(
        "--plot-cell-width",
        action="store_true",
        help="Keep MPAS-Tools cell-width plots in mesh_scratch.",
    )
    parser.add_argument("--force", action="store_true", help="Regenerate existing files.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run_root = args.run_root.expanduser().resolve()
    (run_root / "meshes").mkdir(parents=True, exist_ok=True)
    (run_root / "runs").mkdir(parents=True, exist_ok=True)
    (run_root / "results").mkdir(parents=True, exist_ok=True)

    summaries = [setup_mesh(args, mesh_name) for mesh_name in args.mesh_list]
    manifest = {
        "mesh_type": "uniform_spherical",
        "earth_radius_m": args.earth_radius,
        "ztop_m": args.ztop,
        "meshes": summaries,
    }
    manifest_path = run_root / "mesh_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
