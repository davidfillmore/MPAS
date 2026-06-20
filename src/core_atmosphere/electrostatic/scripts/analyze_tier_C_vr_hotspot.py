#!/usr/bin/env python3
"""Tier C matched-CVT spatial analysis of the FV/TRiSK Poisson solved-phi error.

A variable-resolution JIGSAW mesh is a *general* centroidal-Voronoi mesh (many
non-hexagonal cells), not the icosahedral SCVT (exactly 12 pentagons). So the
"12-pentagon defect" generalises to a topological-defect *field* (nEdgesOnCell!=6),
and the honest control is a *uniform* mesh built on the same fine density grid
(also a general CVT), not the icosahedral A.2 baseline.

This script decomposes the per-cell solved-phi error by
  - region:  core / transition / far  (great-circle distance from the patch centre)
  - defect:  hexagon vs non-hexagon   (nEdgesOnCell == 6 or not)
and reports error-density (sum err^2 / sum area) and local relative-L2 per group,
a within-h-band regression of per-cell error against mesh-distortion metrics, and
the matched VR-vs-uniform-CVT comparison. Reads MPAS output.nc + grid.nc; reuses
the analytic Y_4^2 MMS.

Run (analysis of one mesh):
    ~/miniconda3/envs/mpas/bin/python \
      src/core_atmosphere/electrostatic/scripts/analyze_tier_C_vr_hotspot.py \
      --run-root ~/Data/MPAS/poisson_tier_C_vr --mesh vr_480_120
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import sys

import numpy as np
import netCDF4 as nc

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_tier_A2_sphere_mms as a2run
from setup_tier_C_vr_meshes import great_circle_deg, ALL_MESHES

DEFAULT_PATCH = dict(lat0=0.0, lon0=0.0, r_in_deg=14.0, r_out_deg=34.0)


# ---------------------------------------------------------------------------
# Per-cell error and mesh fields
# ---------------------------------------------------------------------------

def per_cell_error(output_nc):
    """Per-cell volume-weighted err^2 and ref^2 (summed over levels) + geometry.

    Returns dict(err2, ref2, area, h_m, lat_deg, lon_deg). The horizontal-isolated
    source gives level-uniform error, so the level sum is just a constant factor.
    """
    with nc.Dataset(output_nc) as ds:
        phi = a2run.require_level_cell(ds, "phi", "nVertLevels")
        lat = a2run.require_cell(ds, "latCell")
        lon = a2run.require_cell(ds, "lonCell")
        area = a2run.require_cell(ds, "areaCell")
        zgrid = a2run.require_level_cell(ds, "zgrid", "nVertLevelsP1")
    zmid = 0.5 * (zgrid[:-1, :] + zgrid[1:, :])
    dz = zgrid[1:, :] - zgrid[:-1, :]
    z_top = float(np.max(zgrid)) or 1.0
    exact = a2run.phi_exact(lat[None, :], lon[None, :], zmid, z_top)
    vol = area[None, :] * dz
    err2 = np.sum(vol * (phi - exact) ** 2, axis=0)
    ref2 = np.sum(vol * exact ** 2, axis=0)
    return {"err2": err2, "ref2": ref2, "area": np.asarray(area),
            "h_m": np.sqrt(area), "lat_deg": np.degrees(lat), "lon_deg": np.degrees(lon)}


def distortion_metrics(grid_nc):
    """Per-cell mesh-distortion metrics from grid.nc.

    Returns dict(n_edges, nonhex, area_ratio, well_centred):
      nonhex       : nEdgesOnCell != 6 (topological defect)
      area_ratio   : max over cellsOnCell of max(A_i/A_j, A_j/A_i) (cell-size gradient)
      well_centred : min over the cell's edges of dvEdge/dcEdge (TPFA-weight proxy)
    """
    with nc.Dataset(grid_nc) as ds:
        ne = np.asarray(ds.variables["nEdgesOnCell"][:])
        coc = np.asarray(ds.variables["cellsOnCell"][:]) - 1
        eoc = np.asarray(ds.variables["edgesOnCell"][:]) - 1
        area = np.asarray(ds.variables["areaCell"][:])
        dv = np.asarray(ds.variables["dvEdge"][:])
        dc = np.asarray(ds.variables["dcEdge"][:])
    n = area.size
    area_ratio = np.ones(n)
    well_centred = np.full(n, np.inf)
    for i in range(n):
        k = int(ne[i])
        nb = coc[i, :k]
        nb = nb[(nb >= 0) & (nb < n)]
        if nb.size:
            area_ratio[i] = float(np.max(np.maximum(area[i] / area[nb], area[nb] / area[i])))
        ed = eoc[i, :k]
        ed = ed[(ed >= 0) & (ed < dc.size)]
        if ed.size:
            ratio = np.divide(dv[ed], dc[ed], out=np.full(ed.shape, np.inf), where=dc[ed] > 0)
            well_centred[i] = float(np.min(ratio))
    return {"n_edges": ne, "nonhex": ne != 6, "area_ratio": area_ratio, "well_centred": well_centred}


def classify_regions(lat_deg, lon_deg, *, lat0, lon0, r_in_deg, r_out_deg):
    """Return (region, d_deg): region in {core, transition, far} by great-circle d."""
    d = great_circle_deg(lat_deg, lon_deg, lat0, lon0)
    region = np.where(d <= r_in_deg, "core", np.where(d < r_out_deg, "transition", "far"))
    return region, d


# ---------------------------------------------------------------------------
# Aggregations
# ---------------------------------------------------------------------------

def group_error_density(labels, err2, ref2, area):
    """Per-group error-density (sum err^2 / sum area) and local relative-L2.

    Returns list of dict(group, n, err2_frac, area_frac, err_density, rel_l2),
    one per unique label, sorted by label.
    """
    labels = np.asarray(labels)
    tot_e, tot_a = float(err2.sum()), float(area.sum())
    rows = []
    for g in sorted(set(labels.tolist())):
        sel = labels == g
        se2, sr2, sa = float(err2[sel].sum()), float(ref2[sel].sum()), float(area[sel].sum())
        rows.append({"group": g, "n": int(sel.sum()),
                     "err2_frac": se2 / tot_e if tot_e > 0 else 0.0,
                     "area_frac": sa / tot_a if tot_a > 0 else 0.0,
                     "err_density": se2 / sa if sa > 0 else 0.0,
                     "rel_l2": math.sqrt(se2 / sr2) if sr2 > 0 else 0.0})
    return rows


def radial_error_profile(d_deg, err2, ref2, area, *, n_bins):
    """Bin all cells by great-circle distance d; per bin err-density and rel-L2."""
    edges = np.linspace(float(d_deg.min()), float(d_deg.max()) + 1e-9, n_bins + 1)
    idx = np.clip(np.digitize(d_deg, edges) - 1, 0, n_bins - 1)
    out = {k: np.zeros(n_bins) for k in ("d_mid", "err_density", "rel_l2", "count")}
    out["d_mid"] = 0.5 * (edges[:-1] + edges[1:])
    for b in range(n_bins):
        sel = idx == b
        out["count"][b] = int(sel.sum())
        if sel.any():
            se2, sr2, sa = err2[sel].sum(), ref2[sel].sum(), area[sel].sum()
            out["err_density"][b] = se2 / sa if sa > 0 else 0.0
            out["rel_l2"][b] = math.sqrt(se2 / sr2) if sr2 > 0 else 0.0
    return out


def metric_regression(err2, ref2, metric, h_m, *, h_bands):
    """Within fixed local-h bands, Spearman + OLS of per-cell relative error on a
    distortion metric. Controls for resolution so the metric effect is isolated."""
    from scipy.stats import spearmanr
    rel = np.full(err2.shape, np.nan)
    ok = ref2 > 0
    rel[ok] = np.sqrt(err2[ok] / ref2[ok])
    rows = []
    for lo, hi in h_bands:
        sel = ok & (h_m >= lo) & (h_m < hi) & np.isfinite(metric)
        n = int(sel.sum())
        if n >= 8:
            rho = float(spearmanr(metric[sel], rel[sel]).statistic)
            slope = float(np.polyfit(metric[sel], rel[sel], 1)[0])
        else:
            rho, slope = float("nan"), float("nan")
        rows.append({"h_lo_km": lo / 1000.0, "h_hi_km": hi / 1000.0, "n": n,
                     "spearman": rho, "ols_slope": slope})
    return rows


def vr_rate(h_repr, l2):
    """Finest-two and OLS log-log slope; h_repr coarse->fine representative scale."""
    h, e = np.asarray(h_repr, float), np.asarray(l2, float)
    order = np.argsort(-h)
    h, e = h[order], e[order]
    finest_two = math.log(e[-2] / e[-1]) / math.log(h[-2] / h[-1])
    ols = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    return {"finest_two": finest_two, "ols": ols}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def analyze_mesh(run_root, mesh, patch):
    """Load one mesh + its solve; return the per-cell arrays and grouped summaries."""
    run_root = pathlib.Path(run_root).expanduser().resolve()
    out_nc = run_root / "runs" / mesh / "output.nc"
    grid_nc = run_root / "meshes" / mesh / "grid.nc"
    pc = per_cell_error(out_nc)
    met = distortion_metrics(grid_nc)
    region, d = classify_regions(pc["lat_deg"], pc["lon_deg"], **patch)
    defect = np.where(met["nonhex"], "nonhex", "hex")
    rxd = np.array([f"{r}/{x}" for r, x in zip(region, defect)])
    return {"mesh": mesh, "pc": pc, "met": met, "region": region, "d": d,
            "defect": defect,
            "by_region": group_error_density(region, pc["err2"], pc["ref2"], pc["area"]),
            "by_defect": group_error_density(defect, pc["err2"], pc["ref2"], pc["area"]),
            "by_region_defect": group_error_density(rxd, pc["err2"], pc["ref2"], pc["area"])}


def h_bands_from(h_m, n=3):
    """n equal-population local-h bands (meters), as (lo, hi) tuples."""
    qs = np.quantile(h_m, np.linspace(0, 1, n + 1))
    qs[-1] += 1.0
    return [(float(qs[i]), float(qs[i + 1])) for i in range(n)]


def print_summary(res, label):
    print(f"\n=== {label}: {res['mesh']} ===")
    print(" by region      :", " | ".join(
        f"{r['group']}: dens={r['err_density']:.3e} relL2={r['rel_l2']:.3e} n={r['n']} e2%={100*r['err2_frac']:.0f}"
        for r in res["by_region"]))
    print(" by defect      :", " | ".join(
        f"{r['group']}: dens={r['err_density']:.3e} relL2={r['rel_l2']:.3e} n={r['n']} e2%={100*r['err2_frac']:.0f}"
        for r in res["by_defect"]))


def write_csv(res, results_dir):
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / f"tier_C_vr_hotspot_{res['mesh']}.csv"
    pc, met = res["pc"], res["met"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lat_deg", "lon_deg", "d_deg", "h_km", "region", "defect",
                    "n_edges", "area_ratio", "well_centred", "err2", "ref2"])
        for i in range(pc["err2"].size):
            w.writerow([f"{pc['lat_deg'][i]:.4f}", f"{pc['lon_deg'][i]:.4f}",
                        f"{res['d'][i]:.4f}", f"{pc['h_m'][i] / 1000.0:.3f}",
                        res["region"][i], res["defect"][i], int(met["n_edges"][i]),
                        f"{met['area_ratio'][i]:.4f}", f"{met['well_centred'][i]:.4f}",
                        f"{pc['err2'][i]:.6e}", f"{pc['ref2'][i]:.6e}"])
    return path


def global_rel_l2(pc):
    return math.sqrt(float(pc["err2"].sum()) / float(pc["ref2"].sum()))


def rate_mode(run_root, meshes, patch):
    """Solved-phi convergence rate across a nested sequence: global, transition,
    and core relative-L2 vs the representative (fine-core) scale, with slopes."""
    rows = []
    for m in meshes:
        res = analyze_mesh(run_root, m, patch)
        by = {r["group"]: r["rel_l2"] for r in res["by_region"]}
        rows.append({"mesh": m, "h_km": float(ALL_MESHES[m]["fine_km"]),
                     "global": global_rel_l2(res["pc"]),
                     "transition": by.get("transition", float("nan")),
                     "core": by.get("core", float("nan"))})
        r = rows[-1]
        print(f" {m}: h={r['h_km']:g}km  global={r['global']:.3e}  "
              f"transition={r['transition']:.3e}  core={r['core']:.3e}")
    h = [r["h_km"] for r in rows]
    for key in ("global", "transition", "core"):
        vals = [r[key] for r in rows]
        if len(rows) >= 2 and all(v > 0 and math.isfinite(v) for v in vals):
            s = vr_rate(h, vals)
            print(f" slope[{key:10s}] finest_two={s['finest_two']:.3f}  ols={s['ols']:.3f}")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-root", type=pathlib.Path, required=True)
    ap.add_argument("--mesh", default="vr_480_120")
    ap.add_argument("--control", default="u_480", help="Matched uniform-CVT control mesh.")
    ap.add_argument("--patch-json", default=None, help="JSON patch dict; else manifest/default.")
    ap.add_argument("--rate", nargs="+", default=None,
                    help="Nested mesh sequence (coarse->fine) for convergence-rate mode.")
    args = ap.parse_args(argv)

    run_root = args.run_root.expanduser().resolve()
    patch = dict(DEFAULT_PATCH)
    manifest = run_root / "mesh_manifest.json"
    if args.patch_json:
        patch.update(json.loads(args.patch_json))
    elif manifest.exists():
        mp = json.loads(manifest.read_text()).get("patch", {})
        patch.update({k: mp[k] for k in ("lat0", "lon0", "r_in_deg", "r_out_deg") if k in mp})

    if args.rate:
        print(f"=== convergence rate (global / transition / core rel-L2) ===")
        rate_mode(run_root, args.rate, patch)
        return 0

    res = analyze_mesh(run_root, args.mesh, patch)
    print_summary(res, "VR")
    ctrl = analyze_mesh(run_root, args.control, patch)
    print_summary(ctrl, "matched uniform-CVT")

    # Defect cost, isolated on the uniform control (resolution held ~constant):
    cd = {r["group"]: r for r in ctrl["by_defect"]}
    if "hex" in cd and "nonhex" in cd and cd["hex"]["err_density"] > 0:
        print(f"\n[defect cost on {args.control}] err-density nonhex/hex = "
              f"{cd['nonhex']['err_density'] / cd['hex']['err_density']:.2f}x")

    # Matched comparison: VR far-field vs uniform control (both ~coarse general CVT)
    vr_far = next((r for r in res["by_region"] if r["group"] == "far"), None)
    if vr_far:
        print(f"[matched sanity] VR far rel-L2 {vr_far['rel_l2']:.3e} vs "
              f"{args.control} global rel-L2 "
              f"{math.sqrt(ctrl['pc']['err2'].sum() / ctrl['pc']['ref2'].sum()):.3e}")

    # Within-h-band metric regressions on the VR mesh
    h_bands = h_bands_from(res["pc"]["h_m"])
    for name, metric in (("area_ratio", res["met"]["area_ratio"]),
                         ("nonhex", res["met"]["nonhex"].astype(float)),
                         ("inv_well_centred", -res["met"]["well_centred"])):
        rows = metric_regression(res["pc"]["err2"], res["pc"]["ref2"], metric,
                                 res["pc"]["h_m"], h_bands=h_bands)
        print(f"[metric {name:16s}] " + " | ".join(
            f"h~{r['h_lo_km']:.0f}-{r['h_hi_km']:.0f}km rho={r['spearman']:+.2f}(n={r['n']})" for r in rows))

    p1 = write_csv(res, run_root / "results")
    p2 = write_csv(ctrl, run_root / "results")
    print(f"\nWrote {p1}\nWrote {p2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
