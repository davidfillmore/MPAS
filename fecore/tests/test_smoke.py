"""Smoke tests for dolfinx 0.10.0 toolchain de-risking (Task 1).

API adaptations from the brief (dolfinx 0.10.0 specifics):
  1. LinearProblem: `petsc_options_prefix` is a required keyword-only arg.
  2. create_mesh: signature is (comm, cells, e, x) — element before coords.
     The brief had (comm, cells, pts, el), which is reversed.
"""
import numpy as np
import pytest
from mpi4py import MPI


def test_dolfinx_unit_square_poisson_is_second_order():
    """Standard volumetric P1 Poisson MMS on the built-in unit square -> slope ~2."""
    import dolfinx
    import ufl
    from dolfinx.fem.petsc import LinearProblem

    def solve(n):
        msh = dolfinx.mesh.create_unit_square(MPI.COMM_WORLD, n, n)
        V = dolfinx.fem.functionspace(msh, ("Lagrange", 1))
        x = ufl.SpatialCoordinate(msh)
        ue = ufl.sin(ufl.pi * x[0]) * ufl.sin(ufl.pi * x[1])
        f = 2 * ufl.pi**2 * ue
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
        L = f * v * ufl.dx
        facets = dolfinx.mesh.locate_entities_boundary(
            msh, 1, lambda p: np.full(p.shape[1], True)
        )
        dofs = dolfinx.fem.locate_dofs_topological(V, 1, facets)
        bc = dolfinx.fem.dirichletbc(0.0, dofs, V)
        # dolfinx 0.10.0: petsc_options_prefix is required (keyword-only)
        uh = LinearProblem(
            a,
            L,
            bcs=[bc],
            petsc_options_prefix="smoke_volumetric",
            petsc_options={"ksp_type": "cg", "pc_type": "hypre"},
        ).solve()
        err = dolfinx.fem.assemble_scalar(
            dolfinx.fem.form((uh - ue) ** 2 * ufl.dx)
        ) ** 0.5
        return float(err)

    e = [solve(n) for n in (16, 32, 64)]
    slope = np.log(e[-2] / e[-1]) / np.log(2.0)
    assert slope > 1.9, f"unit-square P1 slope {slope:.4f} (expected > 1.9)"


def test_dolfinx_supports_2d_manifold_in_3d():
    """A single triangle embedded in R^3 builds + assembles a stiffness matrix.

    Guards the §7 risk: dolfinx must handle a 2-D mesh with geometric dim 3.

    API adaptation from brief:
      create_mesh signature is (comm, cells, e, x) — element before coords.
    """
    import dolfinx
    import ufl
    import basix

    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 1.0]])  # gdim=3
    cells = np.array([[0, 1, 2]], dtype=np.int64)
    # dolfinx 0.10.0: create_mesh(comm, cells, e, x) — e before x
    el = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
    msh = dolfinx.mesh.create_mesh(MPI.COMM_WORLD, cells, el, pts)
    assert msh.topology.dim == 2 and msh.geometry.dim == 3, (
        f"expected tdim=2, gdim=3; got tdim={msh.topology.dim}, gdim={msh.geometry.dim}"
    )
    V = dolfinx.fem.functionspace(msh, ("Lagrange", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    A = dolfinx.fem.assemble_matrix(
        dolfinx.fem.form(ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx)
    )
    A.scatter_reverse()  # finalize
    assert V.dofmap.index_map.size_global == 3, (
        f"expected 3 global dofs; got {V.dofmap.index_map.size_global}"
    )
