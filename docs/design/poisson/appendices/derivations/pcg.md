(app-pcg)=
# PCG algorithm derivation and convergence

This appendix derives the PCG algorithm from first principles and
records the standard convergence bound.

## Problem setup

Let $A \in \mathbb{R}^{N \times N}$ be symmetric positive
definite and $\bc \in \mathbb{R}^{N}$. We seek
$\bphi \in \mathbb{R}^{N}$ solving

$$
  A \bphi \;=\; \bc.
$$

The $A$-inner product and $A$-norm are

$$
  \langle \vu, \vv \rangle_A \;=\; \vv^{\mathsf{T}} A \vu,
  \qquad
  \|\vu\|_A \;=\; \sqrt{\langle \vu, \vu \rangle_A}.
$$

## Variational formulation

Define

$$
  \mathcal{J}(\vu) \;=\; \tfrac{1}{2} \vu^{\mathsf{T}} A \vu
                         - \vu^{\mathsf{T}} \bc.
$$ (eq-E-J)

The gradient of $\mathcal{J}$ is

$$
  \nabla \mathcal{J}(\vu) \;=\; A \vu - \bc.
$$

Setting $\nabla \mathcal{J} = 0$ yields $A \vu = \bc$: the exact
solution minimizes $\mathcal{J}$. PCG is an algorithm for solving
the minimization problem via successive one-dimensional
minimizations along search directions $\bm{p}^{(k)}$.

## Line search in direction $\bm{p}$

Given current iterate $\bphi$ and direction $\bm{p}$, the next
iterate $\bphi + \alpha \bm{p}$ minimizes $\mathcal{J}$ along the
line when $\alpha$ satisfies

$$
  0 \;=\; \frac{d}{d\alpha}
          \mathcal{J}(\bphi + \alpha \bm{p})
    \;=\; \alpha \bm{p}^{\mathsf{T}} A \bm{p}
          - \bm{p}^{\mathsf{T}} (\bc - A \bphi)
    \;=\; \alpha \bm{p}^{\mathsf{T}} A \bm{p}
          - \bm{p}^{\mathsf{T}} \bm{r},
$$

where $\bm{r} \equiv \bc - A \bphi$ is the residual. Solving for
$\alpha$:

$$
  \alpha \;=\; \frac{\bm{p}^{\mathsf{T}} \bm{r}}
                    {\bm{p}^{\mathsf{T}} A \bm{p}}.
$$ (eq-E-alpha)

## $A$-conjugate search directions

To accelerate convergence, successive search directions are
chosen $A$-conjugate:

$$
  \bm{p}^{(k)\mathsf{T}} A \bm{p}^{(k+1)} \;=\; 0.
$$ (eq-E-conjugate)

With this constraint, after $k + 1$ line searches the iterate
minimizes $\mathcal{J}$ over the Krylov subspace
$\mathcal{K}_{k+1}(A, \bm{r}^{(0)}) =
\mathrm{span}\{\bm{r}^{(0)}, A \bm{r}^{(0)}, \dots, A^{k} \bm{r}^{(0)}\}$,
not only over the most recent direction. This is the key property
that makes CG converge in at most $N$ steps in exact arithmetic.

The next direction is built as a $\beta$-adjustment of the
previous one:

$$
  \bm{p}^{(k+1)} \;=\; \bm{z}^{(k+1)} + \beta_k \bm{p}^{(k)},
$$ (eq-E-pnew)

where $\bm{z}^{(k+1)}$ is the preconditioned residual (see
below) and $\beta_k$ is chosen to enforce
{eq}`eq-E-conjugate`:

$$
  \bm{p}^{(k)\mathsf{T}} A \bm{p}^{(k+1)} \;=\; 0
  \;\;\Rightarrow\;\;
  \beta_k \;=\; -\frac{\bm{p}^{(k)\mathsf{T}} A \bm{z}^{(k+1)}}
                       {\bm{p}^{(k)\mathsf{T}} A \bm{p}^{(k)}}.
$$

Using the residual-update identity
$\bm{r}^{(k+1)} = \bm{r}^{(k)} - \alpha_k A \bm{p}^{(k)}$, a
standard manipulation {cite:p}`shewchuk1994introduction`
rewrites $\beta_k$ in the algorithmic form:

$$
  \beta_k \;=\;
  \frac{\langle \bm{r}^{(k+1)}, \bm{z}^{(k+1)} \rangle}
       {\langle \bm{r}^{(k)},   \bm{z}^{(k)}   \rangle}.
$$

## Preconditioning

Apply $\mathbb{M}^{-1}$ to the residual to get

$$
  \bm{z}^{(k)} \;=\; \mathbb{M}^{-1} \bm{r}^{(k)}.
$$

The preconditioned algorithm operates as if on the transformed
system

$$
  (\mathbb{M}^{-1/2} A \mathbb{M}^{-1/2})
  (\mathbb{M}^{1/2} \bphi) \;=\; \mathbb{M}^{-1/2} \bc.
$$

If $\mathbb{M}^{-1/2} A \mathbb{M}^{-1/2}$ has a smaller
condition number than $A$ alone, PCG converges faster. The art of
preconditioning lies in choosing $\mathbb{M}$ that (i)
approximates $A$ well, (ii) is SPD itself, and (iii) is cheap to
solve.

## Convergence bound

The classical bound for CG applied to SPD $A$ (exact arithmetic)
is {cite:p}`golub2013matrix`

$$
  \|\bphi^{(k)} - \bphi^{*}\|_A
  \;\le\;
  2 \left( \frac{\sqrt{\kappa(A)} - 1}
                 {\sqrt{\kappa(A)} + 1} \right)^{k}
  \|\bphi^{(0)} - \bphi^{*}\|_A,
$$ (eq-E-conv-bound)

where $\kappa(A) = \lambda_{\max}(A) / \lambda_{\min}(A)$ is the
condition number and $\bphi^{*} = A^{-1} \bc$ is the exact
solution. Iteration count to reduce the error by a factor of
$\epsilon$ therefore scales as

$$
  k_\epsilon \;\lesssim\; \tfrac{1}{2}
  \sqrt{\kappa(A)} \log(2/\epsilon).
$$

For the discrete Laplacian on a quasi-uniform mesh of
characteristic spacing $h$, the condition number grows as
$\kappa(A) \sim h^{-2}$, so PCG iteration count grows as
$h^{-1}$. This is the $\mathcal{O}(N^{1/2})$ growth alluded to in
the main-text solver section.

## The three preconditioners implemented

**None ($\mathbb{M} = \mathsf{I}$).** Straight CG. Used as a
baseline reference.

**Jacobi ($\mathbb{M} = \mathsf{D} = \mathrm{diag}(A)$).** Per
iteration, for each $n = (i, k)$:

$$
  z_n \;=\; \mathsf{D}^{-1}_{nn} r_n
       \;=\; \frac{r_n}{|{\mathsf{L}_{nn}}|}.
$$

Cost: one division per cell-column per iteration. Cheap. Reduces
condition number by a modest constant factor but does not change
the $\kappa \sim h^{-2}$ scaling.

**Block symmetric Gauss–Seidel.** Within each MPI rank $p$, let
$\mathcal{I}_p$ denote the rank's owned cells (no halo). Define
the triangular splittings

$$
  A|_{\mathcal{I}_p} \;=\; \mathsf{L}_p + \mathsf{D}_p + \mathsf{U}_p,
$$

with $\mathsf{D}_p$ the rank-$p$ diagonal and
$\mathsf{L}_p, \mathsf{U}_p$ the strictly lower and upper
triangular parts (under the cell-column ordering). The block SGS
preconditioner is

$$
  \mathbb{M}^{-1}
  \;=\; (\mathsf{D}_p + \mathsf{U}_p)^{-1}
        \mathsf{D}_p
        (\mathsf{D}_p + \mathsf{L}_p)^{-1},
$$

applied block-Jacobi across ranks (no coupling between ranks
inside $\mathbb{M}^{-1}$). It is SPD because it is symmetric by
construction and diagonally dominant within each rank. Cost per
iteration: two triangular solves (forward + backward sweep) plus
one halo exchange. Typical iteration-count reduction: 2× on this
problem class.

## Termination criterion

The PCG algorithm terminates when
$\|\bm{r}^{(k)}\|_2 / \|\bc\|_2 < \mathrm{tol}$. Note that
{eq}`eq-E-conv-bound` bounds the $A$-norm of the error, not the
Euclidean norm of the residual; the two are related by

$$
  \|\bm{r}^{(k)}\|_2
  \;=\; \|A(\bphi^{(k)} - \bphi^{*})\|_2
  \;\le\; \|A\|_2 \|\bphi^{(k)} - \bphi^{*}\|_2
  \;\le\; \sqrt{\kappa(A)} \cdot \|A\|_2^{1/2} \|\bphi^{(k)} - \bphi^{*}\|_A,
$$

so a residual tolerance of $10^{-8}$ corresponds to an error
tolerance of roughly $10^{-8} / \sqrt{\kappa(A) \|A\|_2}$, which
for our problems is an excessively tight error bound. In
practice, we use residual convergence as a simple heuristic.
