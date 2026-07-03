#!/usr/bin/env python3
"""Prepare and run a terrain charge-coupled supercell case."""

from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

import netCDF4 as nc
import numpy as np


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_charge_coupled_supercell as charge_supercell  # noqa: E402
import run_terrain_zgrid_mms as terrain_gate  # noqa: E402
import run_tier_A1_cartesian_mms as tier_a1  # noqa: E402
from run_tripole_supercell import seed_run_dir_from_template  # noqa: E402


DEFAULT_RUN_DIR = pathlib.Path(
    "~/Data/MPAS/poisson_charge_coupled_supercell_terrain_h1000/run"
)
DEFAULT_TEMPLATE_RUN_DIR = pathlib.Path("~/Data/MPAS/supercell")
DEFAULT_HILL_HEIGHT = 1000.0
DEFAULT_HILL_HALF_WIDTH_KM = 20.0


def cosine_bell_terrain(xcell, ycell, *, hill_height, half_width_m):
    x0 = 0.5 * (float(np.nanmin(xcell)) + float(np.nanmax(xcell)))
    y0 = 0.5 * (float(np.nanmin(ycell)) + float(np.nanmax(ycell)))
    radius = np.sqrt((xcell - x0) ** 2 + (ycell - y0) ** 2)
    terrain = np.zeros_like(xcell, dtype=float)
    inside = radius < half_width_m
    terrain[inside] = 0.5 * hill_height * (
        1.0 + np.cos(np.pi * radius[inside] / half_width_m)
    )
    return terrain


METRIC_RTOL = 1.0e-10
METRIC_ATOL = 1.0e-9


def _oriented(var, dims_want):
    """Return var's values in canonical dims_want order.

    Init files written by init_atmosphere store (nCells/nEdges, levels)
    C-order; synthetic/legacy files may store the exact reverse. Anything
    else is rejected loudly.
    """
    dims = var.dimensions
    values = np.asarray(var[:])
    if dims == tuple(dims_want):
        return values
    if dims == tuple(reversed(dims_want)):
        return values.transpose(tuple(reversed(range(values.ndim))))
    raise ValueError(f"{var.name} has unsupported dimensions {dims}")


def _write_oriented(var, values, dims_want):
    if var.dimensions == tuple(dims_want):
        var[:] = values
    else:
        var[:] = values.transpose(tuple(reversed(range(values.ndim))))


def read_mesh_arrays(dataset):
    """Horizontal-mesh inputs of the metric formulas (terrain-independent).

    deriv_two comes straight from the file: it is built by
    atm_initialize_advection_rk from horizontal geometry only
    (mpas_init_atm_cases.F:1573), so the flat template's values remain
    valid under any zgrid rewrite.
    """
    v = dataset.variables
    mesh = {
        "cellsOnEdge": np.asarray(
            _oriented(v["cellsOnEdge"], ("nEdges", "TWO")), dtype=np.int64
        ),
        "cellsOnCell": np.asarray(
            _oriented(v["cellsOnCell"], ("nCells", "maxEdges")), dtype=np.int64
        ),
        "nEdgesOnCell": np.asarray(v["nEdgesOnCell"][:], dtype=np.int64),
        "derivTwo": np.asarray(
            _oriented(v["deriv_two"], ("nEdges", "TWO", "FIFTEEN")), dtype=float
        ),
        "dcEdge": np.asarray(v["dcEdge"][:], dtype=float),
        "dvEdge": np.asarray(v["dvEdge"][:], dtype=float),
        "areaCell": np.asarray(v["areaCell"][:], dtype=float),
    }
    n_cells = mesh["areaCell"].size
    if mesh["cellsOnEdge"].min() < 1 or mesh["cellsOnEdge"].max() > n_cells:
        raise ValueError(
            "cellsOnEdge has entries outside 1..nCells; "
            "boundary edges are unsupported (supercell mesh is periodic)"
        )
    return mesh


def recompute_zeta_metrics(dzw):
    """Zeta-only 1-D metrics; mpas_init_atm_cases.F:1629-1662.

    Terrain-independent (zw is the computational coordinate) — recomputed
    only so verify_init_metrics can validate the transcription. Level-1
    entries of dzu/rdzu/fzm/fzp stay 0: the Fortran loops start at k=2 and
    the Registry zero-initializes (template file confirms zeros).
    """
    nz1 = dzw.size
    dzu = np.zeros(nz1)
    rdzu = np.zeros(nz1)
    fzm = np.zeros(nz1)
    fzp = np.zeros(nz1)
    dzu[1:] = 0.5 * (dzw[1:] + dzw[:-1])   # F:1635  dzu(k) = .5*(dzw(k)+dzw(k-1))
    rdzu[1:] = 1.0 / dzu[1:]               # F:1636  rdzu(k) = 1./dzu(k)
    fzp[1:] = 0.5 * dzw[1:] / dzu[1:]      # F:1644  (linear_interpolation)
    fzm[1:] = 0.5 * dzw[:-1] / dzu[1:]     # F:1645
    cof1 = (2.0 * dzu[1] + dzu[2]) / (dzu[1] + dzu[2]) * dzw[0] / dzu[1]  # F:1658
    cof2 = dzu[1] / (dzu[1] + dzu[2]) * dzw[0] / dzu[2]                   # F:1659
    return {
        "dzu": dzu,
        "rdzu": rdzu,
        "fzm": fzm,
        "fzp": fzp,
        "cf1": fzp[1] + cof1,              # F:1660  CF1 = FZP(2) + COF1
        "cf2": fzm[1] - cof1 - cof2,       # F:1661  CF2 = FZM(2) - COF1 - COF2
        "cf3": cof2,                       # F:1662  CF3 = COF2
    }


def recompute_zgrid_metrics(zgrid, dzw, mesh, theta_adv_order):
    """zgrid-derived metrics zz, zxu, zb, zb3 (review finding #9).

    zz, zxu: the supercell/squall-line case, mpas_init_atm_cases.F:1676-1687.
    zb, zb3: the generic terrain formula from the JW case, F:1138-1186 —
    the supercell case writes zeros (F:1962-1963, valid only on a flat
    plane). Only the config_theta_adv_order == 3 branch is implemented
    (the template's stored global attribute).
    """
    if int(theta_adv_order) != 3:
        raise ValueError(
            f"only config_theta_adv_order == 3 is supported, got {theta_adv_order}"
        )
    zgrid = np.asarray(zgrid, dtype=float)
    nz = zgrid.shape[1]
    c1 = mesh["cellsOnEdge"][:, 0] - 1
    c2 = mesh["cellsOnEdge"][:, 1] - 1
    dc = mesh["dcEdge"]
    dv = mesh["dvEdge"]
    area = mesh["areaCell"]
    deriv_two = mesh["derivTwo"]
    coc = mesh["cellsOnCell"]
    nec = mesh["nEdgesOnCell"]

    # zz(k,i) = (zw(k+1)-zw(k)) / (zgrid(k+1,i)-zgrid(k,i))           F:1677
    zz = dzw[None, :] / np.diff(zgrid, axis=1)

    # zxu(k,e) = .5*(zg(k,c2)-zg(k,c1) + zg(k+1,c2)-zg(k+1,c1))/dc    F:1685
    dzg = zgrid[c2, :] - zgrid[c1, :]
    zxu = 0.5 * (dzg[:, :-1] + dzg[:, 1:]) / dc[:, None]

    zb = np.zeros((c1.size, 2, nz))
    zb3 = np.zeros_like(zb)
    for k in range(nz - 1):                # F:1145  do k = 1, nVertLevels
        zk = zgrid[:, k]
        # d2fdx2 = deriv_two(1,side,e)*zg(k,cell) + neighbor terms  F:1153-1165
        d2_1 = deriv_two[:, 0, 0] * zk[c1]
        d2_2 = deriv_two[:, 1, 0] * zk[c2]
        for j in range(coc.shape[1]):
            use1 = (j < nec[c1]) & (coc[c1, j] > 0)
            use2 = (j < nec[c2]) & (coc[c2, j] > 0)
            d2_1[use1] += deriv_two[use1, 0, j + 1] * zk[coc[c1[use1], j] - 1]
            d2_2[use2] += deriv_two[use2, 1, j + 1] * zk[coc[c2[use2], j] - 1]
        z_edge = 0.5 * (zk[c1] + zk[c2]) - dc**2 * (d2_1 + d2_2) / 12.0  # F:1167
        z_edge3 = -(dc**2) * (d2_1 - d2_2) / 12.0                        # F:1171
        zb[:, 0, k] = (z_edge - zk[c1]) * dv / area[c1]                  # F:1178
        zb[:, 1, k] = (z_edge - zk[c2]) * dv / area[c2]                  # F:1179
        zb3[:, 0, k] = z_edge3 * dv / area[c1]                           # F:1180
        zb3[:, 1, k] = z_edge3 * dv / area[c2]                           # F:1181
    # Top-interface row (k = nVertLevels+1) is never written by the
    # Fortran loop and stays 0, matching the stored file.
    return {"zz": zz, "zxu": zxu, "zb": zb, "zb3": zb3}


def verify_init_metrics(init_nc, rtol=METRIC_RTOL, atol=METRIC_ATOL):
    """Golden identity check for the metric transcription (finding #9).

    Recompute every zgrid-derived metric from the file's own zgrid/rdzw
    and compare against the stored fields. On an untouched flat init file
    this validates the Python transcription of the Fortran formulas; on a
    rewritten file it validates self-consistency. dss is excluded: it is
    xnutr*[...] with xnutr = 0. (mpas_init_atm_cases.F:1576,1688-1697),
    identically zero for any zgrid.
    """
    with nc.Dataset(init_nc) as dataset:
        v = dataset.variables
        zgrid = np.asarray(
            _oriented(v["zgrid"], ("nCells", "nVertLevelsP1")), dtype=float
        )
        dzw = 1.0 / np.asarray(v["rdzw"][:], dtype=float).ravel()
        mesh = read_mesh_arrays(dataset)
        theta_adv_order = int(getattr(dataset, "config_theta_adv_order", 3))
        projection = str(
            getattr(dataset, "config_interface_projection", "linear_interpolation")
        ).strip()

        expected = recompute_zeta_metrics(dzw)
        if projection == "layer_integral":  # F:1647-1651 swaps the pair
            expected["fzm"], expected["fzp"] = expected["fzp"], expected["fzm"]
        expected.update(recompute_zgrid_metrics(zgrid, dzw, mesh, theta_adv_order))

        stored = {
            "dzu": np.asarray(v["dzu"][:], dtype=float).ravel(),
            "rdzu": np.asarray(v["rdzu"][:], dtype=float).ravel(),
            "fzm": np.asarray(v["fzm"][:], dtype=float).ravel(),
            "fzp": np.asarray(v["fzp"][:], dtype=float).ravel(),
            "cf1": float(np.asarray(v["cf1"][:])),
            "cf2": float(np.asarray(v["cf2"][:])),
            "cf3": float(np.asarray(v["cf3"][:])),
            "zz": np.asarray(_oriented(v["zz"], ("nCells", "nVertLevels")), dtype=float),
            "zxu": np.asarray(_oriented(v["zxu"], ("nEdges", "nVertLevels")), dtype=float),
            "zb": np.asarray(
                _oriented(v["zb"], ("nEdges", "TWO", "nVertLevelsP1")), dtype=float
            ),
            "zb3": np.asarray(
                _oriented(v["zb3"], ("nEdges", "TWO", "nVertLevelsP1")), dtype=float
            ),
        }

    report = {}
    failures = []
    for name, want in expected.items():
        got = np.asarray(stored[name])
        report[name] = float(np.max(np.abs(got - np.asarray(want))))
        if not np.allclose(got, want, rtol=rtol, atol=atol):
            failures.append(f"{name} (max abs diff {report[name]:.3e})")
    if failures:
        raise RuntimeError(
            f"{init_nc}: stored metrics disagree with recompute: "
            + ", ".join(failures)
        )
    return report


def rewrite_supercell_init_terrain(init_nc, *, hill_height, hill_half_width_m):
    """Rewrite supercell_init.nc with terrain-following columns and metrics.

    Besides the monotone terrain-following zgrid columns and ter, this
    recomputes every zgrid-derived vertical metric the atmosphere core reads
    from the init stream but never re-derives from zgrid — zz, zxu, zb, zb3
    (review finding #9). All writes happen only after everything is computed,
    so a validation failure leaves the file untouched.
    """
    if not np.isfinite(hill_half_width_m) or hill_half_width_m <= 0.0:
        raise ValueError(
            f"terrain half-width must be finite and positive, got {hill_half_width_m}"
        )
    if not np.isfinite(hill_height):
        raise ValueError(f"terrain hill height must be finite, got {hill_height}")

    with nc.Dataset(init_nc, "r+") as dataset:
        xcell = np.asarray(dataset.variables["xCell"][:], dtype=float)
        ycell = np.asarray(dataset.variables["yCell"][:], dtype=float)
        zgrid = dataset.variables["zgrid"]
        dims = zgrid.dimensions
        if "nCells" not in dims or "nVertLevelsP1" not in dims:
            raise ValueError(f"zgrid has unsupported dimensions {dims}")

        cell_axis = dims.index("nCells")
        level_axis = dims.index("nVertLevelsP1")
        if cell_axis == level_axis:
            raise ValueError(f"zgrid has unsupported dimensions {dims}")

        n_cells = zgrid.shape[cell_axis]
        nlevels_p1 = zgrid.shape[level_axis]
        if n_cells != xcell.size or n_cells != ycell.size:
            raise ValueError(
                f"zgrid nCells={n_cells} does not match xCell/yCell size "
                f"{xcell.size}/{ycell.size}"
            )

        ztop = float(np.nanmax(np.asarray(zgrid[:], dtype=float)))
        if ztop <= 0.0:
            raise ValueError(f"zgrid top must be positive, got {ztop}")

        terrain = cosine_bell_terrain(
            xcell,
            ycell,
            hill_height=hill_height,
            half_width_m=hill_half_width_m,
        )
        terrain_min = float(np.nanmin(terrain))
        terrain_max = float(np.nanmax(terrain))
        if (
            not np.isfinite(terrain_min)
            or not np.isfinite(terrain_max)
            or terrain_min < 0.0
            or terrain_max >= ztop
        ):
            raise ValueError(
                "terrain range must be finite and satisfy "
                f"0 <= terrain <= terrain_max < ztop; got "
                f"[{terrain_min}, {terrain_max}] with ztop {ztop}"
            )

        zeta = np.linspace(0.0, ztop, nlevels_p1)
        # Terrain-following column map (run-plots plan Design; identical to
        # the Fortran zgrid with ah(k)=1, hx=ter — mpas_init_atm_cases.F:1673).
        # Shared single source with the MMS gate runner (finding #13).
        znew = terrain_gate.terrain_following_columns(zeta, terrain, ztop)

        # Recompute every zgrid-derived metric (review finding #9): the
        # atmosphere core reads zz/zxu/zb/zb3 from the init stream and
        # never re-derives them from zgrid.
        dzw = 1.0 / np.asarray(dataset.variables["rdzw"][:], dtype=float).ravel()
        if not np.allclose(dzw, ztop / (nlevels_p1 - 1), rtol=1.0e-9):
            raise ValueError(
                "supercell init has a non-uniform zeta grid; the linspace "
                "terrain map assumes the case's hard-coded str = 1.0 "
                "(mpas_init_atm_cases.F:1596)"
            )
        mesh = read_mesh_arrays(dataset)
        theta_adv_order = int(getattr(dataset, "config_theta_adv_order", 3))
        metrics = recompute_zgrid_metrics(znew, dzw, mesh, theta_adv_order)

        _write_oriented(dataset.variables["zgrid"], znew, ("nCells", "nVertLevelsP1"))
        if "ter" in dataset.variables:
            dataset.variables["ter"][:] = terrain
        _write_oriented(dataset.variables["zz"], metrics["zz"], ("nCells", "nVertLevels"))
        _write_oriented(dataset.variables["zxu"], metrics["zxu"], ("nEdges", "nVertLevels"))
        _write_oriented(dataset.variables["zb"], metrics["zb"], ("nEdges", "TWO", "nVertLevelsP1"))
        _write_oriented(dataset.variables["zb3"], metrics["zb3"], ("nEdges", "TWO", "nVertLevelsP1"))

    return {
        "nCells": int(n_cells),
        "ztop": float(ztop),
        "terrain_min": terrain_min,
        "terrain_max": terrain_max,
        "zz_max": float(metrics["zz"].max()),
        "zxu_max_abs": float(np.abs(metrics["zxu"]).max()),
    }


def seed_terrain_run_dir(template_run_dir, run_dir):
    """Seed a terrain run directory, then force supercell_init.nc to be a copy."""
    template_run_dir = template_run_dir.expanduser().resolve()
    run_dir = run_dir.expanduser().resolve()
    if template_run_dir == run_dir:
        raise ValueError(
            f"template_run_dir and run_dir resolve to the same directory: {run_dir}"
        )

    seed_run_dir_from_template(template_run_dir, run_dir)
    source = template_run_dir / "supercell_init.nc"
    target = run_dir / "supercell_init.nc"
    if target.exists() or target.is_symlink():
        target.unlink()
    shutil.copy2(source, target)


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
    hill_height,
    hill_half_width_km,
):
    """Create and configure an isolated terrain dynamic-coupling run directory."""
    seed_terrain_run_dir(template_run_dir, run_dir)
    run_dir = run_dir.expanduser().resolve()
    stats = rewrite_supercell_init_terrain(
        run_dir / "supercell_init.nc",
        hill_height=hill_height,
        hill_half_width_m=hill_half_width_km * 1000.0,
    )
    tier_a1.symlink_model(model.expanduser(), run_dir)

    namelist = run_dir / "namelist.atmosphere"
    streams = run_dir / "streams.atmosphere"
    if not namelist.exists() or not streams.exists():
        raise FileNotFoundError(
            f"{run_dir} needs namelist.atmosphere and streams.atmosphere; "
            f"seed it from a complete supercell run directory"
        )

    charge_supercell.configure_namelist(
        namelist,
        run_duration=run_duration,
        interval=interval,
        stub_alpha=stub_alpha,
        stub_beta=stub_beta,
        poisson_tol=poisson_tol,
        poisson_max_iter=poisson_max_iter,
        solve_at_init=solve_at_init,
    )
    terrain_gate.enable_zgrid_terrain(run_dir / "namelist.atmosphere")
    charge_supercell.configure_streams(streams, output_interval)
    charge_supercell.ensure_coupled_output_stream_list(
        run_dir / "stream_list.atmosphere.output"
    )
    partition = tier_a1.find_partition_file(run_dir, ranks, namelist.read_text())
    if (run_dir / "supercell_init.nc").is_symlink():
        raise RuntimeError(f"{run_dir / 'supercell_init.nc'} must be a real copied file")
    return run_dir, partition, stats


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template-run-dir",
        type=pathlib.Path,
        default=DEFAULT_TEMPLATE_RUN_DIR,
        help="Existing supercell run directory to seed inputs from.",
    )
    parser.add_argument(
        "--run-dir",
        type=pathlib.Path,
        default=DEFAULT_RUN_DIR,
        help="Isolated terrain run directory to create or reuse.",
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
    parser.add_argument("--hill-height", type=float, default=DEFAULT_HILL_HEIGHT)
    parser.add_argument(
        "--hill-half-width-km",
        type=float,
        default=DEFAULT_HILL_HALF_WIDTH_KM,
    )
    parser.add_argument(
        "--verify-init-metrics",
        type=pathlib.Path,
        default=None,
        metavar="INIT_NC",
        help="Recompute zgrid-derived metrics for INIT_NC, compare against "
        "the stored fields, print per-field max diffs, and exit.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.verify_init_metrics is not None:
        report = verify_init_metrics(args.verify_init_metrics.expanduser())
        for name in sorted(report):
            print(f"{name}: max abs diff {report[name]:.3e}")
        print("init metrics verified")
        return 0
    run_dir, partition, stats = prepare_run_dir(
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
        hill_height=args.hill_height,
        hill_half_width_km=args.hill_half_width_km,
    )
    print(f"Prepared {run_dir}")
    print(f"Using partition {partition}")
    print(
        "Terrain range "
        f"{stats['terrain_min']:.3f} to {stats['terrain_max']:.3f} m; "
        f"ztop {stats['ztop']:.3f} m"
    )
    print(
        f"Metrics: zz_max {stats['zz_max']:.6f}, "
        f"|zxu|_max {stats['zxu_max_abs']:.6f}"
    )
    if args.prepare_only:
        return 0
    charge_supercell.run_mpas(run_dir, args.ranks, args.mpiexec)
    print(f"Wrote {run_dir / 'output.nc'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
