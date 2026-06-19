"""Task 6: Surface convergence gate — FEM Poisson on SCVT-dual and icosa meshes.

THE HEADLINE: Does proper FEM (P1-consistent or P2) reach >=1.9 slope on the
SCVT-dual mesh where FV gives ~1.38?  This file defines the gate.

Tests
-----
test_icosa_p2_is_second_order
    icosa family, P2 consistent: slope >= 1.9.

test_scvt_p1_lumped_reproduces_fv_anchor  [skipped if SCVT meshes absent]
    P1-lumped on SCVT-dual should match the FV defect (~1.38); anchors the
    comparison.  Band: 1.2 < slope < 1.6.

test_scvt_proper_fem_beats_the_defect  [skipped if SCVT meshes absent]
    THE HEADLINE: best of P1-consistent and P2-consistent >= 1.9.
    GO  -> leave the assert as-is.
    NO-GO -> xfail with a verdict note; threshold is NEVER lowered.
"""

import pytest
from fecore.verify import convergence


# ---------------------------------------------------------------------------
# icosa gate — always runs
# ---------------------------------------------------------------------------

def test_icosa_p2_is_second_order():
    """P2 FEM on icosa unit-sphere meshes must achieve slope >= 1.9."""
    s = convergence.surface_convergence("icosa", degree=2, rhs_mode="consistent")["slope"]
    assert s >= 1.9, f"icosa P2 slope {s:.4f} < 1.9"


# ---------------------------------------------------------------------------
# SCVT gates — skipped when mesh bundle is absent
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not convergence.scvt_available(),
    reason="SCVT mesh bundle not present at expected path",
)
@pytest.mark.xfail(
    reason=(
        "NO-GO: FEM-P1-lumped on the SCVT-dual achieves slope ~2.00, NOT the "
        "~1.38 FV defect. The ~1.38 defect is a property of the MPAS FV "
        "stencil (finite-volume divergence theorem), not of mass-lumped FEM "
        "on the same mesh. This test's premise — that lumped RHS alone "
        "reproduces the defect — is empirically falsified on these SCVT meshes. "
        "Actual slope = 2.0003."
    ),
    strict=True,
)
def test_scvt_p1_lumped_reproduces_fv_anchor():
    """FEM-P1-lumped on the SCVT-dual should match the FV defect (~1.38).

    Band: 1.2 < slope < 1.6.  This anchors the comparison so we know
    P1-lumped + SCVT really does land in the defect regime.

    VERDICT (NO-GO): Actual slope is ~2.00, not ~1.38.  The FV defect is not
    reproduced by P1-lumped FEM on the same dual mesh.  See report for
    interpretation.
    """
    s = convergence.surface_convergence("scvt", degree=1, rhs_mode="lumped")["slope"]
    assert 1.2 < s < 1.6, f"P1-lumped anchor slope {s:.4f} (expected ~1.38)"


@pytest.mark.skipif(
    not convergence.scvt_available(),
    reason="SCVT mesh bundle not present at expected path",
)
def test_scvt_proper_fem_beats_the_defect():
    """THE HEADLINE: proper FEM (P2 and/or P1-consistent) reaches >=1.9 on the
    SAME SCVT-dual mesh where FV is ~1.38.

    GO  -> this assert passes; do not touch it.
    NO-GO -> mark xfail with a verdict note; NEVER lower 1.9 to make it pass.
    """
    best = max(
        convergence.surface_convergence("scvt", d, m)["slope"]
        for d, m in [(1, "consistent"), (2, "consistent")]
    )
    assert best >= 1.9, (
        f"best proper-FEM SCVT slope {best:.4f} < 1.9 "
        "(proper FEM does NOT recover full 2nd order on SCVT-dual)"
    )
