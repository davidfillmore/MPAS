"""Task 5: Tests for surface Poisson operator (fecore.operator.poisson).

TDD: these tests are written first and must fail until solve_surface is
implemented in fecore/operator/poisson.py.
"""
import numpy as np
import pytest
from fecore.mesh import icosa
from fecore.operator import poisson
from fecore.verify import mms


def test_surface_solve_runs_and_is_spd_consistent():
    """solve_surface returns a finite relative L2 error < 1 for P1-consistent."""
    msh = icosa.icosa_surface_mesh(3)
    uh, err = poisson.solve_surface(
        msh, R=1.0, degree=1, rhs_mode="consistent",
        source_fn=mms.surface_source, exact_fn=mms.y42_cart,
    )
    assert err is not None and 0 < err < 1.0, (
        f"Expected 0 < err < 1.0, got err={err}"
    )


def test_lumped_and_consistent_differ():
    """Consistent vs lumped RHS give different errors (the lever is real)."""
    msh = icosa.icosa_surface_mesh(3)
    _, ec = poisson.solve_surface(msh, 1.0, 1, "consistent",
                                  mms.surface_source, mms.y42_cart)
    _, el = poisson.solve_surface(msh, 1.0, 1, "lumped",
                                  mms.surface_source, mms.y42_cart)
    assert abs(ec - el) > 1e-12, (
        f"Consistent err={ec} and lumped err={el} must differ by > 1e-12"
    )
