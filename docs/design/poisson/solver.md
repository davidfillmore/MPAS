(sec-solver)=
# Solver

## Preconditioned conjugate gradient

The system {eq}`eq-linear-system` is solved using preconditioned
conjugate gradient (PCG) iteration
{cite:p}`hestenes1952methods,saad2003iterative,shewchuk1994introduction`.
The acronym PCG means an iterative conjugate-gradient solver
augmented by a preconditioner that approximates the inverse of the
linear operator. Since $-\mathbb{L}$ is symmetric positive
definite (SPD), the framework solves the equivalent system
$(-\mathbb{L}) \, \boldsymbol{\varphi} = -\boldsymbol{b}$ in
standard form, running PCG with the convention that the operator
applied inside the iteration is $-\mathbb{L}$ and the right-hand
side is $-\boldsymbol{b}$. The pseudocode is written for the SPD
operator $A \equiv -\mathbb{L}$ with right-hand side
$\boldsymbol{c} \equiv -\boldsymbol{b}$:

```text
Given initial guess phi^(0) (default: previous-step phi)
  r^(0) <- c - A phi^(0)
  z^(0) <- M^-1 r^(0)
  p^(0) <- z^(0)
  for k = 0, 1, 2, ...:
      alpha_k    <- <r^(k), z^(k)> / <p^(k), A p^(k)>
      phi^(k+1)  <- phi^(k) + alpha_k p^(k)
      r^(k+1)    <- r^(k) - alpha_k A p^(k)
      if ||r^(k+1)|| / ||c|| < tol: return phi^(k+1), k+1
      z^(k+1)    <- M^-1 r^(k+1)
      beta_k     <- <r^(k+1), z^(k+1)> / <r^(k), z^(k)>
      p^(k+1)    <- z^(k+1) + beta_k p^(k)
```

The tolerance defaults to $10^{-8}$ in relative residual, the
maximum iteration count to $1000$.

## Preconditioners

Three preconditioners are implemented and runtime-selectable:

- **None:** $\mathbb{M} = \mathsf{I}$. Baseline.
- **Jacobi:** $\mathbb{M} = \mathsf{D} \equiv \mathrm{diag}(A)$.
  One division per cell-column per iteration. Default.
- **Block symmetric Gauss–Seidel (block SGS):** within each MPI
  rank's partition, a forward plus a backward SGS sweep;
  Jacobi-style coupling at partition boundaries. Convergence rate
  improves by a constant factor at the cost of one extra halo
  exchange per iteration.

## Null-space projection

Under the Phase 1 boundary conditions the system matrix is
non-singular and the constant-vector null space is absent. A
projection step against the constant vector is nonetheless
implemented behind a runtime flag for use in future pure-Neumann
configurations (e.g. the global current-continuity operator
$\nabla \cdot (\sigma \nabla)$ with conductivity-modulated BCs).

## Halo exchange and parallel implementation

Each matvec $A \, \boldsymbol{p}$ requires one horizontal halo
exchange of $\boldsymbol{p}$; the vertical direction is not
decomposed and requires no halo. MPAS-A's established halo
machinery (`mpas_dmpar_exch_halo_field`) is reused verbatim.

## Warm start

The solution from the previous solve serves as the initial guess
$\boldsymbol{\varphi}^{(0)}$. Because $\rho$ evolves on advection
and sedimentation timescales much longer than a dynamics step,
the warm-started iteration typically converges in a small
fraction of the cold-start iteration count after the initial
transient.
