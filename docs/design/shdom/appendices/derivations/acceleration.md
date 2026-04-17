(app-shdom-accel)=
# Sequence acceleration of the Picard iteration

The Picard iteration {eq}`eq-rte-Picard` converges at a rate
governed by the spectral radius of the composite
scatter-plus-stream operator
$\mathcal{K} = \mathcal{S}\mathcal{L}$. For absorbing
($\omega < 1$) or optically thin media this is comfortably
smaller than unity and the iteration converges in $\sim 10$
passes. For conservative or near-conservative scattering
($\omega \to 1$) in optically thick media — the case of
convective clouds relevant to the SHDOM–MPAS-A coupling — the
spectral radius approaches unity and plain Picard iteration may
require hundreds of passes. This section derives Eqs. (A19)–(A20)
of Evans (1998), which compose each pair of Picard iterates into
an extrapolated estimate designed to cancel the leading
geometric-convergence error.

(app-shdom-accel-model)=
## Iteration difference vector and 2D geometric model

Let $J^{(n)}$ denote the spherical-harmonic source-function field
at iteration $n$, regarded as a vector in the
$(N_{\mathrm{gridpoints}} \times N_{\mathrm{lm}})$-dimensional
space of source-function states; and define the *iteration
difference* vector

$$
  u^{(n)} \;\equiv\; J^{(n)} - J^{(n-1)}.
$$ (eq-u-diff)

Under the assumption that the iteration converges geometrically
with a *single* dominant (possibly complex) rate, the difference
vectors spiral toward zero:

$$
  u^{(n)} \;=\; r\,\hat{R}_\psi\,u^{(n-1)},
$$ (eq-u-geometric)

with $r = |u^{(n)}|/|u^{(n-1)}|$ the per-iteration length-
contraction ratio and $\hat{R}_\psi$ a rotation by angle $\psi$
in the 2D plane spanned by $u^{(n-1)}$ and $u^{(n)}$. Both $r$
and $\psi$ are measured directly from the two most recent
iterates,

$$
  \boxed{\;
    r \;=\; \frac{|u^{(n)}|}{|u^{(n-1)}|},
    \qquad
    \cos\psi
    \;=\; \frac{u^{(n-1)}\cdot u^{(n)}}{|u^{(n-1)}|\,|u^{(n)}|},
  \;}
$$ (eq-A19)

which is Eq. (A19) of Evans (1998). The 2D reduction is
justified by the observation that the iteration lives in the
two-dimensional subspace spanned by the two most recent
difference vectors, within which the dominant geometric mode can
be resolved; finer spectral structure of $\mathcal{K}$ is
neglected, consistent with the convergence-acceleration
literature this method draws on.

(app-shdom-accel-A20)=
## Extrapolation distance: Eq. (A20)

Under the geometric model {eq}`eq-u-geometric`, the infinite sum
of remaining differences from iteration $n - 1$ is

$$
  \sum_{k=0}^{\infty} u^{(n-1+k)}
  \;=\; \sum_{k=0}^{\infty} r^{k}\,\hat{R}_{\psi}^{k}\,u^{(n-1)}
  \;=\; \bigl(\mathsf{I} - r\,\hat{R}_{\psi}\bigr)^{-1}\,u^{(n-1)}.
$$ (eq-geom-sum-2D)

In the 2D plane where the spiral lives,
$\mathsf{I} - r\hat{R}_\psi$ has determinant
$1 + r^{2} - 2r\cos\psi$ and its inverse is a scalar multiple of
a rotated copy of the identity; the *length* of the summed
vector is

$$
  \bigl|\bigl(\mathsf{I} - r\hat{R}_\psi\bigr)^{-1}\,u^{(n-1)}\bigr|
  \;=\; \frac{|u^{(n-1)}|}{\sqrt{1 + r^{2} - 2r\cos\psi}}.
$$ (eq-sum-length)

Evans (1998) extrapolates the source function by a scalar
multiple $a$ of the *current* difference $u^{(n)}$:

$$
  J^{\prime\,(n)}_{\mathrm{accel}}
  \;=\; J^{(n)} + a\,\bigl(J^{(n)} - J^{(n-1)}\bigr)
  \;=\; J^{(n)} + a\,u^{(n)},
$$ (eq-accel-update)

with $a$ chosen so that the extrapolated point lies at the
geometric-series limit along the $u^{(n-1)}$ direction. Evans's
formula for $a$ is

$$
  \boxed{\;
    a \;=\; \frac{1 - r\cos\psi + r^{\,1 + \pi/(2\psi)}}
                 {1 + r^{2} - 2r\cos\psi}
          \;-\; 1,
  \;}
$$ (eq-A20)

reproducing Eq. (A20) of Evans (1998). Two structural
observations on {eq}`eq-A20` guide its interpretation:

- The denominator is exactly the scalar
  $\det(\mathsf{I} - r\hat{R}_\psi) = 1 + r^{2} - 2r\cos\psi$
  encountered in {eq}`eq-geom-sum-2D`. Its role is to normalize
  the extrapolation by the geometric-series contraction factor,
  so $a \to 0$ as $r \to 0$ (no acceleration needed for
  fast-converging iterations) and $a$ diverges only when $r = 1$
  and $\psi = 0$ (the identity-operator limit where no progress
  is made at all).
- The numerator contains the unusual term
  $r^{1 + \pi/(2\psi)}$, which encodes the "pitch" of the 2D
  geometric spiral: a spiral with rotation angle $\psi$ completes
  a half-turn in $\pi/\psi$ iterations, and
  $r^{\pi/(2\psi)}$ is the length contraction over a
  quarter-turn. This heuristic factor, not present in a naive
  geometric-series sum {eq}`eq-sum-length`, is the 2D analogue
  of the Aitken-$\Delta^{2}$ correction applied to a spiraling
  rather than a purely contracting sequence.

We do not derive {eq}`eq-A20` from first principles here; Evans
(1998) presents it as a heuristic adopted from the
convergence-acceleration literature cited in his appendix. In
practice, $a$ from {eq}`eq-A20` yields a factor of 2–4 reduction
in iteration count for conservative-scattering cases
(Evans, 1998, appendix §e), and the acceleration is applied
*every other* iteration to allow the geometric-convergence
estimate to stabilize between applications. Because the
accelerated update can overshoot the true fixed point in
pathological cases, SHDOM permits a configuration-time toggle to
disable sequence acceleration entirely.
