(app-implementation-note)=
# Implementation note: volume-integrated vs per-unit-volume form

This appendix compares two mathematically equivalent forms of the
discrete operator and documents the one chosen for
implementation.

## Two equivalent forms

**Form I: volume-integrated.**
$\mathbb{L}_{\mathrm{I}} \bphi = \bb_{\mathrm{I}}$ with

$$
\begin{aligned}
  (\mathbb{L}_{\mathrm{I}} \bphi)_{i, k}
  &\;=\; \sum_e \varepsilon \frac{\elen \dz{k}}{\dcen} (\varphi_j - \varphi_i)
         + \varepsilon A_i \left[
           \tfrac{\varphi_{k+1} - \varphi_k}{\dz{k+\hlf}}
           + \tfrac{\varphi_{k-1} - \varphi_k}{\dz{k-\hlf}}
           \right], \\
  b_{\mathrm{I}, i, k} &\;=\; -\rho_{i, k} \vol.
\end{aligned}
$$ (eq-H-form1)

$\mathbb{L}_{\mathrm{I}}$ is symmetric in standard $\ell^{2}$
inner product (see [](symmetry.md)). CG uses standard inner
products and norms:
$\langle \vu, \vv \rangle = \sum_n u_n v_n$.

**Form II: per-unit-volume.**
$\mathbb{L}_{\mathrm{II}} \bphi = \bb_{\mathrm{II}}$ with

$$
\begin{aligned}
  (\mathbb{L}_{\mathrm{II}} \bphi)_{i, k}
  &\;=\; \frac{1}{\vol} (\mathbb{L}_{\mathrm{I}} \bphi)_{i, k}, \\
  b_{\mathrm{II}, i, k} &\;=\; -\rho_{i, k}.
\end{aligned}
$$

$\mathbb{L}_{\mathrm{II}}$ is *not* symmetric in standard
$\ell^{2}$ because the diagonal scaling by $\vol^{-1}$ breaks
symmetry at the matrix level. It is, however, symmetric in the
volume-weighted inner product
$\langle \vu, \vv \rangle_V = \sum_n V_n u_n v_n$:

$$
\begin{aligned}
  \langle \mathbb{L}_{\mathrm{II}} \vu, \vv \rangle_V
  &\;=\; \sum_n V_n \cdot (\mathbb{L}_{\mathrm{I}} \vu)_n / V_n \cdot v_n \\
  &\;=\; \sum_n (\mathbb{L}_{\mathrm{I}} \vu)_n v_n \\
  &\;=\; \langle \mathbb{L}_{\mathrm{I}} \vu, \vv \rangle_{\ell^{2}},
\end{aligned}
$$

and similarly
$\langle \vu, \mathbb{L}_{\mathrm{II}} \vv \rangle_V =
\langle \mathbb{L}_{\mathrm{I}} \vv, \vu \rangle_{\ell^{2}} =
\langle \mathbb{L}_{\mathrm{I}} \vu, \vv \rangle_{\ell^{2}}$
(since $\mathbb{L}_{\mathrm{I}}$ is $\ell^{2}$-symmetric). CG in
Form II uses $V$-weighted inner products and norms throughout.

**Equivalence.** Both forms solve the same continuous PDE to the
same order. The two systems have identical solutions $\bphi$. The
choice of form affects how the matrix and the inner products are
expressed in code but not the mathematical result.

## This narrative uses Form I

The main text and previous appendices use Form I throughout,
because:

1. Form I's matrix is $\ell^{2}$-symmetric, so the generic CG
   algorithm derived in [](pcg.md) applies with no modification.
2. The symmetry proof in [](symmetry.md) is direct and
   constructive.
3. Form I's right-hand side is $-\rho \vol$, which is physically
   interpretable as total charge in the cell.

Form II is a valid alternative but introduces the bookkeeping of
a $V$-weighted inner product and is less transparent
pedagogically.

## Which form the implementation uses

The implementation adopts **Form I**. With Form I, all CG code
operates in standard $\ell^{2}$ inner products (no $V$-weighted
variants needed), the matrix is transparently symmetric, and the
right-hand side $b = -\rho \cdot \vol$ has the physical
interpretation of total charge in the cell-column.

Concretely, Form I fixes the stored operator weights as:

- `h_weight(k, e)` $= \varepsilon \cdot \elen \cdot \dz{k} / \dcen$
  — a 2D array over `(nVertLevels, nEdges)` because $\dz{k}$
  depends on $k$ (and, on terrain-following meshes, on $i$ as
  well).
- `v_weight_upper(k, i)` $= \varepsilon \cdot A_i / \dz{k+\hlf}$,
  zero at $k = K$ under the Neumann top condition.
- `v_weight_lower(k, i)` $= \varepsilon \cdot A_i / \dz{k-\hlf}$
  for $k > 1$ and $= \varepsilon \cdot A_i / (\hlf \dz{1})$ at
  $k = 1$ to represent the ground-Dirichlet face.

The matvec for row $(i, k)$ assembles the sum of face-flux
contributions directly from these weights and returns
$A \, \varphi$ where $A \equiv -\mathbb{L}$ is the SPD negated
Laplacian used inside the CG iteration.
