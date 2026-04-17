(app-symmetry)=
# Symmetry and positive-definiteness of $\mathbb{L}$

This appendix proves Property (P2) of the main text:
$-\mathbb{L}$ is symmetric positive definite.

## Matrix entries

Let us index the unknowns by a flat index $n = (i, k)$ ranging
over all $N = n_{\mathrm{cells}} \times n_{\mathrm{levels}}$
cell-columns. Denote the matrix $\mathbb{L}$ by $\mathsf{L}$ with
entries $\mathsf{L}_{nm}$. From [](horizontal-operator.md) and
[](vertical-operator.md), the non-zero entries are:

$$
\begin{aligned}
  \mathsf{L}_{(i, k), (j, k)}
    &\;=\; +\varepsilon_e \, \frac{\elen \dz{k}}{\dcen}
    \qquad \text{for } j = \mathrm{nbr}(i, e), \\
  \mathsf{L}_{(i, k), (i, k+1)}
    &\;=\; +\varepsilon \, A_i / \dz{k+\hlf}, \\
  \mathsf{L}_{(i, k), (i, k-1)}
    &\;=\; +\varepsilon \, A_i / \dz{k-\hlf}, \\
  \mathsf{L}_{(i, k), (i, k)}
    &\;=\; -\sum_{e \in \partial i} \varepsilon_e \elen \dz{k}/\dcen
           - \varepsilon A_i / \dz{k+\hlf}
           - \varepsilon A_i / \dz{k-\hlf}
           + \Delta^{\mathrm{bc}}_{i, k}.
\end{aligned}
$$ (eq-C-L-entries)

In the diagonal entry, $\Delta^{\mathrm{bc}}_{i, k}$ carries the
boundary correction:
$\Delta^{\mathrm{bc}}_{i, 1} = -\varepsilon A_i / \dz{\tfrac{1}{2}}^{\mathrm{half}} + \varepsilon A_i / \dz{\tfrac{1}{2}}$
(the second term cancels the interior $\dz{k-\hlf}$ contribution
since there is no layer $k = 0$); for $k = K$,
$\Delta^{\mathrm{bc}}_{i, K} = +\varepsilon A_i / \dz{K+\hlf}$
cancels the interior upper-face term. For interior $k$,
$\Delta^{\mathrm{bc}}_{i, k} = 0$.

## Symmetry: $\mathsf{L}_{nm} = \mathsf{L}_{mn}$

We show $\mathsf{L}_{nm} = \mathsf{L}_{mn}$ for every pair
$(n, m)$.

**Horizontal pair.** For $n = (i, k)$, $m = (j, k)$ with
$j = \mathrm{nbr}(i, e)$:

$$
  \mathsf{L}_{(i, k), (j, k)} \;=\; \varepsilon_e \elen \dz{k} / \dcen.
$$

By the consistency result of
{ref}`app-horizontal-consistency`, $\elen$, $\dcen$, $\dz{k}$,
and $\varepsilon_e$ are all intrinsic to the shared edge $e$ and
do not depend on whether we approach the edge from cell $i$ or
cell $j$. Hence

$$
  \mathsf{L}_{(j, k), (i, k)} \;=\; \varepsilon_e \elen \dz{k} / \dcen
  \;=\; \mathsf{L}_{(i, k), (j, k)}.
$$

**Vertical pair.** For $n = (i, k)$, $m = (i, k+1)$:

$$
  \mathsf{L}_{(i, k), (i, k+1)} \;=\; \varepsilon A_i / \dz{k+\hlf}.
$$

The "upper face of layer $k$" is geometrically identical to the
"lower face of layer $k + 1$"; both have area $A_i$, and the
midpoint-to-midpoint distance
$\dz{k+\hlf} = z_{i, k+1} - z_{i, k}$ is the same. Therefore

$$
  \mathsf{L}_{(i, k+1), (i, k)} \;=\; \varepsilon A_i / \dz{k+\hlf}
  \;=\; \mathsf{L}_{(i, k), (i, k+1)}.
$$

**Non-adjacent pair.** If cells $n$ and $m$ share no face, both
$\mathsf{L}_{nm}$ and $\mathsf{L}_{mn}$ are zero.

Since all non-zero entries come in symmetric pairs and diagonal
entries are trivially equal to themselves, $\mathsf{L}$ is
symmetric.

## Negative-semidefiniteness for interior rows

The diagonal entry of row $n = (i, k)$ is minus the sum of the
absolute values of off-diagonal entries in that row, plus a
boundary correction:

$$
  \mathsf{L}_{nn} \;=\; -\sum_{m \ne n} |\mathsf{L}_{nm}|
                        + \Delta^{\mathrm{bc}}_{i, k}.
$$

For interior rows, $\Delta^{\mathrm{bc}}_{i, k} = 0$, so the row
has zero row sum: $\sum_m \mathsf{L}_{nm} = 0$. This means the
constant vector $\bm{1}$ is a null vector of $\mathsf{L}$ on the
interior rows.

By the Gershgorin circle theorem, each eigenvalue $\lambda$ of
$\mathsf{L}$ lies in the union of disks centred at the diagonal
entries $\mathsf{L}_{nn}$ with radii
$\sum_{m \ne n} |\mathsf{L}_{nm}|$. For zero-row-sum rows with
negative diagonal, each Gershgorin disk is contained in
$\{\lambda : \lambda \le 0\}$, so all eigenvalues satisfy
$\lambda \le 0$.

## Strict negative-definiteness under ground Dirichlet

With the ground-Dirichlet condition, the $k = 1$ rows acquire
$\Delta^{\mathrm{bc}}_{i, 1} = -\varepsilon A_i /
\dz{\tfrac{1}{2}}^{\mathrm{half}} + \varepsilon A_i /
\dz{\tfrac{1}{2}} \ne 0$. Specifically, the second term cancels
the interior $k - \hlf$ contribution (which for $k = 1$ would
have pointed at the nonexistent $k = 0$), and the first term is
the negative contribution $-\varepsilon A_i /
\dz{\tfrac{1}{2}}^{\mathrm{half}}$ from the Dirichlet face.
Combined, the $k = 1$ diagonal picks up
$-\varepsilon A_i / \dz{\tfrac{1}{2}}^{\mathrm{half}}$ relative
to a row whose diagonal would have summed to zero.

Therefore
$\sum_m \mathsf{L}_{(i, 1), m} = -\varepsilon A_i / \dz{\tfrac{1}{2}}^{\mathrm{half}} < 0$:
the $k = 1$ rows are strictly diagonally dominant. The remaining
rows have zero row sum.

Equivalently, the Gershgorin disks for $k = 1$ rows are strictly
contained in $\{\lambda < 0\}$. This rules out $\lambda = 0$ as
an eigenvalue, because an eigenvector for $\lambda = 0$ would
satisfy $\mathsf{L} \bm{v} = 0$ everywhere, including on $k = 1$
rows; but a strictly diagonally dominant row cannot annihilate a
non-zero vector whose $(i, 1)$-entries differ from those of its
neighbours. A short argument: suppose $\mathsf{L} \bm{v} = 0$.
On $k = 1$:

$$
  -\left( \text{sum of off-diags} + \varepsilon A_i /\dz{\tfrac{1}{2}}^{\mathrm{half}} \right) v_{i, 1}
  + \sum_{m \sim (i, 1)} |\mathsf{L}_{(i, 1), m}|\, v_m \;=\; 0.
$$

If $v$ attains its maximum value at some $(i, 1)$, then
$v_m \le v_{i, 1}$ for all neighbours $m$, and the above sum
yields $v_{i, 1} \cdot (-\varepsilon A_i / \dz{\tfrac{1}{2}}^{\mathrm{half}}) \ge 0$,
which requires $v_{i, 1} \le 0$. By symmetry (if $v$ attains
minimum at $(i', 1)$), $v_{i', 1} \ge 0$. So $v$ at $k = 1$ is
bounded by $[0, 0]$; i.e. $v$ vanishes on $k = 1$. By a
discrete-maximum-principle argument applied recursively upward in
$k$, $v$ vanishes everywhere, contradicting non-triviality.

Hence $\lambda = 0$ is not an eigenvalue, and combined with
$\lambda \le 0$ for all eigenvalues from the Gershgorin argument,
$\mathsf{L}$ has strictly negative eigenvalues. Equivalently,
$-\mathsf{L}$ is symmetric positive definite.

## Implication for PCG

A standard result of Krylov iteration
{cite:p}`golub2013matrix,saad2003iterative` is that conjugate
gradients converge for any SPD operator. The PCG algorithm
(see [](pcg.md)) applied with $A = -\mathsf{L}$ and
$\bm{c} = -\bb$ therefore converges to the unique solution
$\bphi$ of $\mathsf{L} \bphi = \bb$.
