"""Surface and shell Poisson operators: FEM Laplace(-Beltrami) solves.

Provides two solvers:

solve_surface — Laplace-Beltrami Poisson on a closed 2-sphere surface mesh.
    Assembles and solves eps * a(u,v) = L(v) on a tdim=2, gdim=3 mesh.
    Gauge: one DOF pinned; error norms use zero-mean alignment.

solve_shell — Full 3-D Laplace Poisson on a spherical shell mesh.
    Assembles and solves eps * a(u,v) = L(v) on a tdim=3, gdim=3 shell mesh.
    BCs: ground Dirichlet (phi=0) on the inner sphere (r≈R); top Neumann
    (zero-flux) on the outer sphere (r≈R+H) — natural, no term needed;
    laterally closed (no lateral boundary on a spherical shell).
    No gauge fix needed: the Dirichlet BC removes the null space, making
    the system SPD.

Both solvers support two RHS modes:
  consistent -- standard assembled load integral(rho v dx).
  lumped      -- mass-lumped point x volume: row-sums of the consistent mass
                 matrix applied to nodal source values; mimics the FV load
                 and reveals the pentagon defect on SCVT-dual meshes.

dolfinx 0.10.0 API notes (confirmed in Task 1 + Task 3 + Task 8):
  - functionspace(mesh, ("Lagrange", degree))
  - locate_entities_boundary(mesh, dim, marker_fn) — marker_fn receives (3,n)
  - locate_dofs_topological(V, entity_dim, entities)
  - dirichletbc(value, dofs, V)
  - assemble_scalar(form(...))
  - LinearProblem(a, L, bcs=[...], petsc_options_prefix=..., petsc_options=...)
  - interpolate: dolfinx passes coords as shape (3, n); transpose to (n, 3)
    before calling user functions with signature f(xyz, R, eps).
"""

import numpy as np
from mpi4py import MPI
import dolfinx
import dolfinx.fem
import dolfinx.fem.petsc
import dolfinx.mesh as _dmesh
import ufl
import petsc4py.PETSc as PETSc


def solve_surface(mesh, R, degree, rhs_mode, source_fn, exact_fn=None, eps=1.0):
    """Solve the Laplace-Beltrami Poisson on a closed sphere-surface mesh.

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh
        Triangulated sphere surface (tdim=2, gdim=3); no boundary.
    R : float
        Sphere radius (passed through to source_fn and exact_fn).
    degree : {1, 2}
        Lagrange polynomial degree.
    rhs_mode : {"consistent", "lumped"}
        How to assemble the right-hand side.
        "consistent": standard integral(rho * v dx).
        "lumped": row-summed mass matrix applied to nodal source values.
    source_fn : callable
        Source (RHS forcing) function with signature
        ``source_fn(xyz, R, eps) -> ndarray``
        where ``xyz`` has shape (n, 3).
    exact_fn : callable or None
        Exact solution, ``exact_fn(xyz, R) -> ndarray``.  Used to pin DOF 0
        and compute the relative L2 error.  If None, DOF 0 is pinned to 0
        and (uh, None) is returned.
    eps : float
        Diffusion coefficient (uniform scalar).  Default 1.0.

    Returns
    -------
    uh : dolfinx.fem.Function
        Discrete solution (zero-mean-aligned when exact_fn is provided).
    l2_rel : float or None
        Relative L2 error ``||uh - u_exact||_0 / ||u_exact||_0`` after
        area-weighted mean subtraction, or None if exact_fn is None.
    """
    V = dolfinx.fem.functionspace(mesh, ("Lagrange", degree))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    # ------------------------------------------------------------------ #
    # Interpolate source into a Function (dolfinx passes x as shape (3,n))
    # ------------------------------------------------------------------ #
    src = dolfinx.fem.Function(V)
    src.interpolate(lambda x: source_fn(x.T, R, eps))

    # ------------------------------------------------------------------ #
    # Exact solution and gauge: pin DOF 0
    # ------------------------------------------------------------------ #
    u_exact = None
    if exact_fn is not None:
        u_exact = dolfinx.fem.Function(V)
        u_exact.interpolate(lambda x: exact_fn(x.T, R))
        pin_value = float(u_exact.x.array[0])
    else:
        pin_value = 0.0

    dofs_bc = np.array([0], dtype=np.int32)
    bc = dolfinx.fem.dirichletbc(np.float64(pin_value), dofs_bc, V)

    # ------------------------------------------------------------------ #
    # Bilinear form (shared between both modes)
    # ------------------------------------------------------------------ #
    a_ufl = eps * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    # ------------------------------------------------------------------ #
    # Assemble and solve
    # ------------------------------------------------------------------ #
    if rhs_mode == "consistent":
        L_ufl = src * v * ufl.dx
        from dolfinx.fem.petsc import LinearProblem
        prob = LinearProblem(
            a_ufl, L_ufl, bcs=[bc],
            petsc_options_prefix="fecore_surf",
            petsc_options={"ksp_type": "cg", "pc_type": "hypre"},
        )
        uh = prob.solve()

    elif rhs_mode == "lumped":
        uh = _solve_lumped(mesh, V, u, v, a_ufl, src, bc, eps)

    else:
        raise ValueError(f"rhs_mode must be 'consistent' or 'lumped', got {rhs_mode!r}")

    # ------------------------------------------------------------------ #
    # Zero-mean alignment + relative L2 error
    # ------------------------------------------------------------------ #
    if exact_fn is None:
        return uh, None

    one_fn = dolfinx.fem.Constant(mesh, dolfinx.default_scalar_type(1.0))
    area = dolfinx.fem.assemble_scalar(dolfinx.fem.form(one_fn * ufl.dx))

    mean_uh = dolfinx.fem.assemble_scalar(dolfinx.fem.form(uh * ufl.dx)) / area
    mean_ex = dolfinx.fem.assemble_scalar(dolfinx.fem.form(u_exact * ufl.dx)) / area

    uh.x.array[:] -= mean_uh
    u_exact.x.array[:] -= mean_ex

    diff2 = dolfinx.fem.assemble_scalar(
        dolfinx.fem.form((uh - u_exact) ** 2 * ufl.dx)
    )
    norm2 = dolfinx.fem.assemble_scalar(
        dolfinx.fem.form(u_exact ** 2 * ufl.dx)
    )
    l2_rel = float(np.sqrt(diff2 / norm2))
    return uh, l2_rel


def solve_shell(mesh, R, H, degree, rhs_mode, source_fn, exact_fn, eps=1.0):
    """Solve the full 3-D Poisson equation on a spherical shell mesh.

    Assembles and solves

        -eps * div(grad(u)) = f   in V (shell)

    subject to:
        u = 0   on the inner sphere r ≈ R  (ground Dirichlet, homogeneous)
        n·grad(u) = 0   on the outer sphere r ≈ R+H  (top Neumann, natural)

    No gauge fix is needed: the Dirichlet BC pins the solution, making the
    system SPD (unlike the closed-surface pure-Neumann case in solve_surface).

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh
        Tetrahedral spherical shell mesh (tdim=3, gdim=3).
    R : float
        Inner sphere radius.
    H : float
        Shell thickness (outer radius = R + H).
    degree : {1, 2}
        Lagrange polynomial degree.
    rhs_mode : {"consistent", "lumped"}
        How to assemble the right-hand side.
        "consistent": standard integral(rho * v dx).
        "lumped": row-summed mass matrix applied to nodal source values.
    source_fn : callable
        Forcing function with signature ``source_fn(xyz, R, H, eps) -> ndarray``
        where ``xyz`` has shape (n, 3).
    exact_fn : callable
        Exact solution with signature ``exact_fn(xyz, R, H) -> ndarray``
        where ``xyz`` has shape (n, 3).
    eps : float
        Diffusion coefficient (uniform scalar).  Default 1.0.

    Returns
    -------
    uh : dolfinx.fem.Function
        Discrete solution on the shell.
    l2_rel : float
        Relative L2 error ||uh - u_exact||_0 / ||u_exact||_0.
    """
    V = dolfinx.fem.functionspace(mesh, ("Lagrange", degree))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    # ------------------------------------------------------------------ #
    # Ground Dirichlet BC: phi = 0 on the inner surface (r ≈ R)
    # Locate inner-surface facets by centroid radius < R + tol.
    # Use 5% of R as tolerance to absorb chord shortfall of flat tets.
    # ------------------------------------------------------------------ #
    tol = 0.05 * R
    inner_facets = _dmesh.locate_entities_boundary(
        mesh, mesh.topology.dim - 1,
        lambda x: np.linalg.norm(x, axis=0) < R + tol,
    )
    inner_dofs = dolfinx.fem.locate_dofs_topological(
        V, mesh.topology.dim - 1, inner_facets
    )
    bc = dolfinx.fem.dirichletbc(np.float64(0.0), inner_dofs, V)

    # ------------------------------------------------------------------ #
    # Interpolate source.  dolfinx passes x as shape (3, n); transpose.
    # ------------------------------------------------------------------ #
    src = dolfinx.fem.Function(V)
    src.interpolate(lambda x: source_fn(x.T, R, H, eps))

    # ------------------------------------------------------------------ #
    # Exact solution (for error computation).
    # ------------------------------------------------------------------ #
    u_exact = dolfinx.fem.Function(V)
    u_exact.interpolate(lambda x: exact_fn(x.T, R, H))

    # ------------------------------------------------------------------ #
    # Bilinear form: full 3-D gradient (not Laplace-Beltrami).
    # ------------------------------------------------------------------ #
    a_ufl = eps * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    # ------------------------------------------------------------------ #
    # Assemble and solve
    # ------------------------------------------------------------------ #
    if rhs_mode == "consistent":
        L_ufl = src * v * ufl.dx
        from dolfinx.fem.petsc import LinearProblem
        prob = LinearProblem(
            a_ufl, L_ufl, bcs=[bc],
            petsc_options_prefix="fecore_shell",
            petsc_options={"ksp_type": "cg", "pc_type": "hypre"},
        )
        uh = prob.solve()

    elif rhs_mode == "lumped":
        uh = _solve_shell_lumped(mesh, V, u, v, a_ufl, src, bc, eps)

    else:
        raise ValueError(f"rhs_mode must be 'consistent' or 'lumped', got {rhs_mode!r}")

    # ------------------------------------------------------------------ #
    # Relative L2 error (no zero-mean alignment — Dirichlet pins gauge)
    # ------------------------------------------------------------------ #
    diff2 = dolfinx.fem.assemble_scalar(
        dolfinx.fem.form((uh - u_exact) ** 2 * ufl.dx)
    )
    norm2 = dolfinx.fem.assemble_scalar(
        dolfinx.fem.form(u_exact ** 2 * ufl.dx)
    )
    l2_rel = float(np.sqrt(diff2 / norm2))
    return uh, l2_rel


def _solve_shell_lumped(mesh, V, u, v, a_ufl, src, bc, eps):
    """Solve the shell system with lumped-mass RHS at the PETSc level.

    1. Assemble stiffness A = eps * integral( grad(u).grad(v) dx ), apply BC.
    2. Assemble consistent mass M = integral( u*v dx ).
    3. Compute lumped mass mL = M * 1  (row sums).
    4. Build load b_i = mL_i * rho_i  (nodal source).
    5. Apply lifting and BC to b; CG-solve A uh = b.
    """
    a_form = dolfinx.fem.form(a_ufl)

    # Stiffness with BC diagonal
    A = dolfinx.fem.petsc.assemble_matrix(a_form, bcs=[bc])
    A.assemble()

    # Consistent mass (no BCs — we need the full lumped volumes)
    m_form = dolfinx.fem.form(u * v * ufl.dx)
    M = dolfinx.fem.petsc.assemble_matrix(m_form)
    M.assemble()

    # Lumped mass: row sums of M
    ones = M.createVecRight()
    ones.set(1.0)
    mL_vec = M.createVecRight()
    M.mult(ones, mL_vec)
    mL = mL_vec.getArray().copy()

    # Load vector b_i = mL_i * rho_i
    rho_nodal = src.x.array.copy()
    b_arr = mL * rho_nodal

    b = M.createVecRight()
    b.setArray(b_arr)
    b.assemble()

    # Apply BC: lifting then set
    dolfinx.fem.petsc.apply_lifting(b, [a_form], bcs=[[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    dolfinx.fem.petsc.set_bc(b, [bc])

    # KSP solve
    ksp = PETSc.KSP().create(mesh.comm)
    ksp.setOperators(A)
    ksp.setType("cg")
    ksp.getPC().setType("hypre")
    ksp.setFromOptions()
    ksp.setUp()

    x_vec = A.createVecRight()
    ksp.solve(b, x_vec)

    if not ksp.is_converged:
        reason = ksp.getConvergedReason()
        raise RuntimeError(
            f"KSP did not converge: reason={reason}. "
            "Check that the system is well-posed (inner Dirichlet applied)."
        )

    uh = dolfinx.fem.Function(V)
    uh.x.array[:] = x_vec.getArray()
    uh.x.scatter_forward()

    # Cleanup PETSc objects
    ksp.destroy()
    A.destroy()
    M.destroy()
    b.destroy()
    x_vec.destroy()
    ones.destroy()
    mL_vec.destroy()

    return uh


def _solve_lumped(mesh, V, u, v, a_ufl, src, bc, eps):
    """Solve with lumped-mass RHS at the PETSc matrix/vector level.

    1. Assemble stiffness A = eps * integral( grad(u).grad(v) dx ), apply BC.
    2. Assemble consistent mass M = integral( u*v dx ).
    3. Compute lumped mass mL = M * 1  (row sums).
    4. Build load b_i = mL_i * rho_i  (nodal source).
    5. Apply lifting and BC to b; CG-solve A uh = b.
    """
    a_form = dolfinx.fem.form(a_ufl)

    # Stiffness with BC diagonal
    A = dolfinx.fem.petsc.assemble_matrix(a_form, bcs=[bc])
    A.assemble()

    # Consistent mass (no BCs — we need the full lumped areas)
    m_form = dolfinx.fem.form(u * v * ufl.dx)
    M = dolfinx.fem.petsc.assemble_matrix(m_form)
    M.assemble()

    # Lumped mass: row sums of M
    ones = M.createVecRight()
    ones.set(1.0)
    mL_vec = M.createVecRight()
    M.mult(ones, mL_vec)
    mL = mL_vec.getArray().copy()

    # Load vector b_i = mL_i * rho_i
    rho_nodal = src.x.array.copy()
    b_arr = mL * rho_nodal

    b = M.createVecRight()
    b.setArray(b_arr)
    b.assemble()

    # Apply BC: lifting then set
    dolfinx.fem.petsc.apply_lifting(b, [a_form], bcs=[[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    dolfinx.fem.petsc.set_bc(b, [bc])

    # KSP solve
    ksp = PETSc.KSP().create(mesh.comm)
    ksp.setOperators(A)
    ksp.setType("cg")
    ksp.getPC().setType("hypre")
    ksp.setFromOptions()
    ksp.setUp()

    x_vec = A.createVecRight()
    ksp.solve(b, x_vec)

    if not ksp.is_converged:
        reason = ksp.getConvergedReason()
        raise RuntimeError(
            f"KSP did not converge: reason={reason}. "
            "Check that the system is well-posed (single DOF pinned)."
        )

    uh = dolfinx.fem.Function(V)
    uh.x.array[:] = x_vec.getArray()
    uh.x.scatter_forward()

    # Cleanup PETSc objects
    ksp.destroy()
    A.destroy()
    M.destroy()
    b.destroy()
    x_vec.destroy()
    ones.destroy()
    mL_vec.destroy()

    return uh
