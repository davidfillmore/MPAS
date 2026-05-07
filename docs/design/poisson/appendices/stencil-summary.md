(app-stencil-summary)=
# Stencil coefficients summary

For each cell-column $(i, k)$, the row of $\mathbb{L}$ from
{eq}`eq-L-full` has the non-zero entries shown below.

| Column | Coefficient |
|---|---|
| $(j, k)$, $j = \mathrm{nbr}(i, e)$ | $\varepsilon \, \ell_e \Delta z_k / d_e$ |
| $(i, k+1)$ | $\varepsilon \, A_i / \Delta z_{k+1/2}$ |
| $(i, k-1)$ | $\varepsilon \, A_i / \Delta z_{k-1/2}$ |
| $(i, k)$ (diagonal) | $-\sum_{e \in \partial i} \varepsilon \, \ell_e \Delta z_k / d_e \,-\, \varepsilon A_i / \Delta z_{k+1/2} \,-\, \varepsilon A_i / \Delta z_{k-1/2}$ |

**Boundary adjustments.** For $k = 1$ rows, the
$\Delta z_{k-1/2}$ term is replaced by the half-layer distance to
the ground face, and the resulting off-diagonal contribution is
absorbed into the diagonal (Dirichlet $\varphi = 0$ at the face
contributes nothing to row $i$'s other entries). For $k = K$
rows, the $\Delta z_{k+1/2}$ term is dropped entirely (Neumann).
