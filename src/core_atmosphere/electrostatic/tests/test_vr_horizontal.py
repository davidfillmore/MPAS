#!/usr/bin/env python3
"""Tests for Task B4: variable-resolution / mesh-distortion characterization.

Compares baseline diagonal-Hodge vs cotangent operator on perturbed
(non-well-centred) planar Voronoi meshes using SOLUTION ERROR (not
operator-residual) as the common yardstick.

Run with:
    ~/miniconda3/envs/mpas/bin/python -m unittest \
        src.core_atmosphere.electrostatic.tests.test_vr_horizontal -v
"""

import importlib.util
import pathlib
import unittest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    REPO_ROOT
    / "src/core_atmosphere/electrostatic/scripts/prototype_vr_horizontal.py"
)


def load():
    spec = importlib.util.spec_from_file_location("prototype_vr_horizontal", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class VRHorizontalTests(unittest.TestCase):

    def test_unperturbed_both_operators_second_order(self):
        """Unperturbed (perturb=0) both operators give solution-error slope > 1.8."""
        m = load()
        slope_base = m.vr_solution_error(0.0, "baseline")
        slope_cot = m.vr_solution_error(0.0, "cotangent")
        print(f"\n[VR unperturbed] baseline={slope_base:.3f} cotangent={slope_cot:.3f}")
        self.assertGreater(slope_base, 1.8)
        self.assertGreater(slope_cot, 1.8)

    def test_cotangent_spd_at_each_perturbation(self):
        """Cotangent operator is SPD (lambda_min > 0) at every perturbation level."""
        m = load()
        for pf in (0.0, 0.1, 0.25, 0.5):
            mesh = m.perturbed_mesh(16, pf, seed=0)
            A = m.vr_operator(mesh, "cotangent")
            n = A.shape[0]
            # Ground by removing row/col 0 (same grounding used in solution-error solve)
            gr = A[1:n, 1:n].tocsr()
            lam_min = m.spd_min_eig(gr)
            print(f"\n[VR SPD] perturb={pf:.2f} lambda_min={lam_min:.4e}")
            self.assertGreater(lam_min, 0.0, f"cotangent not SPD at perturb={pf}")

    def test_sweep_records_both_operators(self):
        """Sweep records baseline and cotangent solution-error slopes for all perturbations."""
        m = load()
        rows = m.vr_sweep()
        self.assertEqual(len(rows), 4)
        for r in rows:
            self.assertIn("baseline_slope", r)
            self.assertIn("cotangent_slope", r)
            self.assertIn("perturb", r)
            print(
                f"\n[VR] perturb={r['perturb']:.2f} "
                f"baseline={r['baseline_slope']:.3f} "
                f"cotangent={r['cotangent_slope']:.3f}"
            )


if __name__ == "__main__":
    unittest.main()
