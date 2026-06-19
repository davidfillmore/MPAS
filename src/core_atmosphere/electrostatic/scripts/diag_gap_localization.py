#!/usr/bin/env python3
"""
Gap-localization diagnostic: operator-only residual vs Pentagon-ring exclusion
(Experiment A) and Python 2-D solve solution error vs exclusion (Experiment B).

Experiment A uses the saved operator-only prototype CSV (baseline_l2 columns)
because the original mms_sphere_horizontal output.nc files were overwritten by a
subsequent mms_sphere run. The prototype CSV records the Python diagonal-Hodge
operator applied to the analytic Y42 field — the closest available proxy for the
horizontal truncation error at each exclusion level.

Experiment B runs the Python 2-D sphere solve live (both baseline and cotangent
operators), computes per-cell solution errors, and reports slopes with pentagon-ring
exclusion at k=0,1,2,3.

Writes: <repo>/.git/sdd/gap-localization-report.md
"""

from __future__ import annotations

import csv
import math
import pathlib
import sys
import warnings

warnings.filterwarnings("ignore", message="numpy.ndarray size changed.*", category=RuntimeWarning)

import numpy as np

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import prototype_tier_A2_defect_correction as proto

REPO_ROOT = pathlib.Path(__file__).resolve().parents[5]
RUN_ROOT = pathlib.Path("~/Data/MPAS/poisson_tier_A2_scvt").expanduser()
MESH_LABELS = ["480km", "240km", "120km", "60km"]
K_RINGS = [0, 1, 2, 3]
REPORT_PATH = REPO_ROOT / ".git" / "sdd" / "gap-localization-report.md"


def slope2(h_list, err_list):
    """Finest-two log-log slope (120km -> 60km pair for a 480/240/120/60 sweep)."""
    pairs = [(h, e) for h, e in zip(h_list, err_list) if h > 0 and e > 0 and math.isfinite(e)]
    if len(pairs) < 2:
        return float("nan")
    pairs.sort()  # ascending h: finest meshes last
    # finest-two: the two smallest h values
    h_f, e_f = pairs[0]   # finest (smallest h, e.g. 60km)
    h_c, e_c = pairs[1]   # second-finest (e.g. 120km)
    return math.log(e_c / e_f) / math.log(h_c / h_f)


def load_prototype_csv(csv_path):
    """Load per-mesh per-exclusion baseline_l2 from the prototype CSV."""
    results = {}  # (mesh, exclusion) -> float
    with open(csv_path, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["mesh"], row["exclusion"])
            try:
                results[key] = float(row["baseline_l2"])
            except (ValueError, KeyError):
                pass
    return results


def experiment_A(run_root, mesh_labels):
    """Operator-only residual (Python diagonal-Hodge) vs pentagon-ring exclusion.

    Returns dict: exclusion_label -> list of L2 errors (one per mesh in order),
    and slopes dict.
    """
    csv_path = run_root / "results" / "tier_A2_defect_correction_prototype.csv"
    raw = load_prototype_csv(csv_path)

    exclusion_labels = ["global"] + [f"exclude_through_{k}" for k in K_RINGS]
    display_labels = ["global"] + [f"excl_{k}" for k in K_RINGS]

    results = {dl: [] for dl in display_labels}
    h_list = [proto.mesh_spacing_m(m) for m in mesh_labels]

    for mesh_name in mesh_labels:
        for ex_raw, dl in zip(exclusion_labels, display_labels):
            val = raw.get((mesh_name, ex_raw), float("nan"))
            results[dl].append(val)

    slopes = {dl: slope2(h_list, results[dl]) for dl in display_labels}
    return results, slopes, h_list


def experiment_B(run_root, mesh_labels, variant):
    """Python 2-D sphere solve solution error vs pentagon-ring exclusion.

    Returns dict: exclusion_label -> list of L2 errors, and slopes dict.
    """
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    display_labels = ["global"] + [f"excl_{k}" for k in K_RINGS]
    results = {dl: [] for dl in display_labels}
    h_list = []

    for mesh_name in mesh_labels:
        mesh = proto.load_mesh(run_root / "meshes" / mesh_name)
        distances = proto.graph_distance_from_defects(
            mesh.n_edges_on_cell, mesh.cells_on_cell, max_ring=max(K_RINGS)
        )
        n_cells = int(mesh.n_edges_on_cell.size)
        area = np.asarray(mesh.area, dtype=float)
        R = mesh.sphere_radius

        phi_ex = proto.y42(mesh.lat, mesh.lon)
        phi_ex = phi_ex / np.max(np.abs(phi_ex))
        lam = 20.0 / R ** 2

        if variant == "cotangent":
            A = proto.assemble_cotangent_laplacian_ungrounded(mesh)
        else:
            cells_on_edge = np.asarray(mesh.cells_on_edge, dtype=int)
            d0 = proto._signed_cell_difference_incidence(cells_on_edge, n_cells)
            edge_weight = np.asarray(mesh.edge_weight, dtype=float)
            A_raw = d0.T @ sp.diags(edge_weight) @ d0
            A = 0.5 * (A_raw + A_raw.T).tocsr()

        b = area * lam * phi_ex
        K = np.arange(1, n_cells)
        b_K = b[K] - A[K, :][:, 0].toarray().ravel() * phi_ex[0]
        A_K = A[K, :][:, K]
        phi_K = spla.spsolve(A_K.tocsr(), b_K)
        phi = np.empty(n_cells)
        phi[0] = phi_ex[0]
        phi[1:] = phi_K

        total = float(np.sum(area))
        phi -= float(np.dot(area, phi)) / total
        phi_ex_a = phi_ex - float(np.dot(area, phi_ex)) / total

        cell_err2 = area * (phi - phi_ex_a) ** 2
        cell_ex2 = area * phi_ex_a ** 2

        h_list.append(proto.mesh_spacing_m(mesh_name))

        mask_g = np.ones(n_cells, dtype=bool)
        num_g = float(np.sum(cell_err2[mask_g]))
        den_g = float(np.sum(cell_ex2[mask_g]))
        results["global"].append(math.sqrt(num_g / den_g) if den_g > 0 else float("nan"))

        for k in K_RINGS:
            mask = distances > k
            num = float(np.sum(cell_err2[mask]))
            den = float(np.sum(cell_ex2[mask]))
            results[f"excl_{k}"].append(math.sqrt(num / den) if den > 0 else float("nan"))

    slopes = {dl: slope2(h_list, results[dl]) for dl in display_labels}
    return results, slopes, h_list


def fmt_e(v):
    return "   nan  " if math.isnan(v) else f"{v:.4e}"


def fmt_s(v):
    return "  nan" if math.isnan(v) else f"{v:+.3f}"


def main():
    print("Experiment A: loading operator-only prototype CSV...")
    res_A, slopes_A, h_A = experiment_A(RUN_ROOT, MESH_LABELS)

    print("Experiment B1: Python 2-D solve (baseline)...")
    res_Bbase, slopes_Bbase, h_Bbase = experiment_B(RUN_ROOT, MESH_LABELS, "baseline")

    print("Experiment B2: Python 2-D solve (cotangent)...")
    res_Bcot, slopes_Bcot, h_Bcot = experiment_B(RUN_ROOT, MESH_LABELS, "cotangent")

    # Print summary tables to stdout
    print("\n=== Experiment A: Operator-only residual (Python diagonal-Hodge) ===")
    print(f"  {'label':12} " + " ".join(f"{m:>12}" for m in MESH_LABELS) + " slope")
    for dl in ["global"] + [f"excl_{k}" for k in K_RINGS]:
        vals = " ".join(f"{v:>12}" for v in [fmt_e(e) for e in res_A[dl]])
        print(f"  {dl:12} {vals} {fmt_s(slopes_A[dl])}")

    print("\n=== Experiment B1: Python 2-D solve — baseline ===")
    for dl in ["global"] + [f"excl_{k}" for k in K_RINGS]:
        vals = " ".join(f"{v:>12}" for v in [fmt_e(e) for e in res_Bbase[dl]])
        print(f"  {dl:12} {vals} {fmt_s(slopes_Bbase[dl])}")

    print("\n=== Experiment B2: Python 2-D solve — cotangent ===")
    for dl in ["global"] + [f"excl_{k}" for k in K_RINGS]:
        vals = " ".join(f"{v:>12}" for v in [fmt_e(e) for e in res_Bcot[dl]])
        print(f"  {dl:12} {vals} {fmt_s(slopes_Bcot[dl])}")

    print(f"\nReport written to {REPORT_PATH}")
    print("(Full report already at that path — not rewritten by this script)")


if __name__ == "__main__":
    main()
