(app-shdom-cellint)=
# Integration along a discrete ordinate across a grid cell

Step 2 of SHDOM's Picard cycle
(see {ref}`app-shdom-rte-picard`) integrates the radiative
transfer equation in its integral form {eq}`eq-eq2` along every
discrete ordinate through every cell of the computational grid.
The integration is performed *one cell at a time*, sweeping the
ordinate direction from the upstream cell face (where the
entering radiance $I(0)$ is known by bilinear interpolation of the
four face-grid values that the ray intersects) to the grid point
at the opposite side of the cell. This section derives the
cell-integration formulas Eqs. (A12)–(A17) of Evans (1998): the
linear-extinction and linear-$(Jk)$ parameterizations within a
cell, the constant-extinction exact integral, and Evans's
small-$\tau$ and hybrid approximations for the general
linear-extinction case.

(app-shdom-cellint-setup)=
## Cell-integration setup and notation

Let $s \in [0, s_{\mathrm{cell}}]$ be arclength along the discrete
ordinate, with $s = 0$ at the upstream cell face (entering the
cell) and $s = s_{\mathrm{cell}}$ at the grid point at the
opposite face. To reduce clutter we write $s$ for
$s_{\mathrm{cell}}$ and let subscript $0$ denote values at the
entering face and subscript $1$ denote values at the grid point
(not to be confused with the angular indices $l, m$ of
[](spherical-harmonics.md) and [](scattering.md)):

$$
  k_0 \equiv k(s' = 0), \quad k_1 \equiv k(s' = s),
  \qquad J_0 \equiv J(s' = 0), \quad J_1 \equiv J(s' = s).
$$

Gridpoint values of $(k, J)$ are supplied by trilinear
interpolation of the property and source fields; face values are
supplied by bilinear interpolation across the two-dimensional
cell face pierced by the ray. Within the cell, SHDOM parameterizes
$k$ and the source function product $(Jk)$ as linearly varying
along the ray.

(app-shdom-cellint-linear-k)=
## Linear extinction model and optical path: Eqs. (A12)–(A13)

Assume

$$
  \boxed{\;
    k(s') \;=\; k_0 + (k_1 - k_0)\,\frac{s'}{s},
  \;}
$$ (eq-A12)

reproducing Eq. (A12) of Evans (1998). The full-cell optical
depth is

$$
  \tau \;\equiv\; \tau(0, s) \;=\; \int_{0}^{s} k(s')\,ds'
              \;=\; \frac{(k_0 + k_1)\,s}{2},
$$ (eq-A13)

which is Eq. (A13), the trapezoidal-rule integral of a linear
function. The partial optical path from an interior point $s'$ to
the grid point $s$ follows by direct integration of {eq}`eq-A12`:

$$
  \tau(s', s) \;=\; \int_{s'}^{s} k(s'')\,ds''
              \;=\; (s - s')\,\Bigl[\frac{k_0 + k_1}{2}
                 \;+\;\frac{(k_1 - k_0)}{2s}\,s'\Bigr].
$$ (eq-tau-partial)

At the endpoints $\tau(0, s) = \tau$ and $\tau(s, s) = 0$ as
required, and differentiation reproduces $-d\tau(s', s)/ds' =
k(s')$.

(app-shdom-cellint-linear-Jk)=
## Linear source–extinction product: Eq. (A14)

SHDOM treats the *product* $(Jk)$, rather than $J$ alone, as
linear across the cell:

$$
  \boxed{\;
    (Jk)(s') \;=\; J_0 k_0 + (J_1 k_1 - J_0 k_0)\,\frac{s'}{s},
  \;}
$$ (eq-A14)

which is Eq. (A14). The rationale was noted at the end of
{ref}`app-shdom-rte-int`: $(Jk)$ is the physical *emission
coefficient per unit path length* (per the $k$-inside-the-integrand
structure of Eq. {eq}`eq-eq2`), and is the natural
linear-interpolation variable. Treating $J$ alone as linear would
leave $(Jk) = J(s') k(s')$ as a *quadratic* in $s'$, complicating
the integrals below.

(app-shdom-cellint-A15)=
## Constant-extinction exact integral: Eq. (A15)

The special case $k_0 = k_1 \equiv k$ admits a closed-form
integral. Under this assumption $\tau(s', s) = k(s - s') = \tau -
k s'$, and the cell-integration formula {eq}`eq-eq2` reads

$$
  I(s) \;=\; e^{-\tau}\,I(0)
  \;+\; \int_{0}^{s} e^{-k(s - s')}\, k\,J(s')\,ds',
$$ (eq-A15-step0)

where $J(s') = J_0 + (J_1 - J_0)\,s'/s$ is linear (equivalent to
{eq}`eq-A14` at constant $k$). Change variables $u \equiv k(s - s')$,
so $du = -k\,ds'$ and $s' = s - u/k = s(1 - u/\tau)$, with bounds
$u\colon \tau \to 0$ as $s'\colon 0 \to s$. Then

$$
\begin{aligned}
  \int_{0}^{s} e^{-k(s - s')}\,k\,J(s')\,ds'
  &\;=\; \int_{0}^{\tau} e^{-u}\,J\bigl(s(1 - u/\tau)\bigr)\,du \\[4pt]
  &\;=\; \int_{0}^{\tau} e^{-u}\Bigl[J_1 + (J_0 - J_1)\,\frac{u}{\tau}\Bigr]\,du,
\end{aligned}
$$

where we used $J(s(1 - u/\tau)) = J_0 + (J_1 - J_0)(1 - u/\tau) =
J_1 + (J_0 - J_1)\,u/\tau$. The two remaining elementary integrals
are

$$
\begin{aligned}
  \int_{0}^{\tau} e^{-u}\,du       &\;=\; 1 - e^{-\tau}, \\
  \int_{0}^{\tau} u\,e^{-u}\,du &\;=\; \bigl[-u\,e^{-u}\bigr]_{0}^{\tau}
                                      + \int_{0}^{\tau} e^{-u}\,du
                                    \;=\; 1 - (1 + \tau)\,e^{-\tau}.
\end{aligned}
$$

Substituting,

$$
\begin{aligned}
  \int_{0}^{\tau} e^{-u}\Bigl[J_1 + (J_0 - J_1)\tfrac{u}{\tau}\Bigr]\,du
  &\;=\; J_1\,(1 - e^{-\tau})
     \;+\; \frac{J_0 - J_1}{\tau}\,\bigl[1 - (1 + \tau)\,e^{-\tau}\bigr] \\[4pt]
  &\;=\; (1 - e^{-\tau})\,J_1
     \;+\; (J_0 - J_1)\,\Bigl[\frac{1 - e^{-\tau}}{\tau} - e^{-\tau}\Bigr],
\end{aligned}
$$

where the second line rearranges
$[1 - (1 + \tau) e^{-\tau}]/\tau = (1 - e^{-\tau})/\tau -
e^{-\tau}$. Combining with the transmitted boundary term,

$$
  \boxed{\;
    I(s) \;=\; e^{-\tau}\,I(0) \;+\; (1 - e^{-\tau})\,J_1
             \;+\; (J_0 - J_1)\,\Bigl[\frac{1 - e^{-\tau}}{\tau} - e^{-\tau}\Bigr],
  \;}
$$ (eq-A15)

reproducing Eq. (A15) of Evans (1998). Three limiting sanity
checks:

- **Constant source** ($J_0 = J_1 \equiv J$): the third term
  vanishes and {eq}`eq-A15` reduces to
  $I(s) = e^{-\tau} I(0) + (1 - e^{-\tau})\,J$, the textbook
  Schwarzschild form for a homogeneous slab.
- **Optically thin** ($\tau \to 0$):
  $(1 - e^{-\tau})/\tau \to 1$ and $e^{-\tau} \to 1$, so the
  bracket $(1 - e^{-\tau})/\tau - e^{-\tau} \to 0$ and
  $(1 - e^{-\tau}) J_1 \to \tau J_1$; expanding to next order
  gives $\tau(J_0 + J_1)/2 - \tau^{2}(2 J_0 + J_1)/6 +
  O(\tau^3)$, the trapezoidal-plus-correction small-$\tau$ limit.
- **Optically thick** ($\tau \to \infty$): $e^{-\tau} \to 0$ and
  $(1 - e^{-\tau})/\tau \to 0$, so $I(s) \to J_1$. The emergent
  radiance saturates at the local grid-point source function, as
  expected when the cell is opaque and the boundary contribution
  is fully absorbed.

(app-shdom-cellint-taylor)=
## Leading Taylor expansion for linear extinction

With linear $k$ and linear $(Jk)$, the full integral in
{eq}`eq-eq2` is not elementary: the exponential of
{eq}`eq-tau-partial` contains a quadratic in $s'$. For
orientation — and to anchor the approximation schemes of the next
subsection — we record the exact leading Taylor expansion in
small cell optical depth.

Change variables $x \equiv s'/s \in [0, 1]$. Then

$$
  \tau(s', s) \;=\; \frac{s}{2}\,\bigl[(1 - x)^{2}\,k_0
                      + (1 - x^{2})\,k_1\bigr],
$$ (eq-tau-dimensionless)

and

$$
  \int_{0}^{s} e^{-\tau(s', s)}\,(Jk)(s')\,ds'
  \;=\; s\,\int_{0}^{1} e^{-\tau(s', s)}\,\bigl[q_0 + \Delta q\,x\bigr]\,dx,
$$

with $q_0 \equiv J_0 k_0$, $q_1 \equiv J_1 k_1$, and
$\Delta q \equiv q_1 - q_0$. For small cell optical depth
$\tau \sim s \max(k_0, k_1) \ll 1$ the exponential expands as
$e^{-\tau(s', s)} = 1 - \tau(s', s) +
\tfrac{1}{2}\tau(s', s)^{2} - \cdots$. Using the elementary
integrals

$$
  \int_{0}^{1}(1 - x)^{2}\,dx = \tfrac{1}{3}, \quad
  \int_{0}^{1}(1 - x)^{2}\,x\,dx = \tfrac{1}{12}, \quad
  \int_{0}^{1}(1 - x^{2})\,dx = \tfrac{2}{3}, \quad
  \int_{0}^{1}(1 - x^{2})\,x\,dx = \tfrac{1}{4},
$$

the two leading orders evaluate to

$$
\begin{aligned}
  O(s): \quad
  s\int_{0}^{1}(q_0 + \Delta q\,x)\,dx
  &\;=\; s\,(q_0 + \tfrac{1}{2}\Delta q)
  \;=\; \tfrac{s}{2}\,(J_0 k_0 + J_1 k_1), \\[4pt]
  O(s^{2}): \quad
  -\,s\int_{0}^{1}\tau(s', s)(q_0 + \Delta q\,x)\,dx
  &\;=\; -\,\tfrac{s^{2}}{24}\bigl[3 k_0 q_0 + 5 k_1 q_0 + k_0 q_1 + 3 k_1 q_1\bigr].
\end{aligned}
$$

Combined with the transmitted-boundary term expanded to the same
order,

$$
  I(s) \;=\; \bigl[1 - \tau + \tfrac{1}{2}\tau^{2}\bigr]\,I(0)
       \;+\; \frac{s}{2}(J_0 k_0 + J_1 k_1)
       \;-\; \frac{s^{2}}{24}\bigl[3 k_0 q_0 + 5 k_1 q_0
              + k_0 q_1 + 3 k_1 q_1\bigr]
       \;+\; O(s^{3}).
$$ (eq-I-taylor)

In the constant-extinction limit $k_0 = k_1 = k$ (so $\tau = k s$),
equation {eq}`eq-I-taylor` reduces to $I(s) = [1 - \tau +
\tau^{2}/2]\,I(0) + \tau(J_0 + J_1)/2 - \tau^{2}(2 J_0 + J_1)/6 +
O(\tau^{3})$, matching the $\tau \to 0$ expansion of {eq}`eq-A15`
noted in sanity check (ii) above. This cross-check confirms the
algebra.

(app-shdom-cellint-A16A17)=
## Evans's closed-form approximations: Eqs. (A16)–(A17)

For linear extinction with $k_0 \ne k_1$, the exact integral is
non-elementary. Evans (1998) uses two closed-form approximations
engineered for computational efficiency:

- Eq. (A16), a compact "expansion solution" designed to be
  accurate at small-to-moderate cell optical depth
  ($\tau \lesssim 2$) and to be computable from a handful of
  arithmetic operations per ray segment.
- Eq. (A17), a hybrid formula that combines (A16) with the
  constant-extinction exact form {eq}`eq-A15`, chosen to recover
  the correct $\tau \to \infty$ asymptote $I(s) \to J_1$ and used
  in SHDOM for all cells with $\tau > 2$.

Evans states (Eq. (A16) of Evans, 1998):

$$
  I(s) \;=\; e^{-\tau}\,I(0) + (1 - e^{-\tau})\,\mathcal{A}_{16},
$$ (eq-A16-schematic)

where $\mathcal{A}_{16}$ is an explicit algebraic combination of
$J_0 k_0$, $J_1 k_1$, $k_0$, $k_1$, and $s$, built to reduce to
a trapezoidal source-function average in the small-$\tau$ limit
and to remain stable in the linear-$k$ case. The hybrid Eq. (A17)
has the same outer structure but with $\mathcal{A}_{17}$
constructed to match $\mathcal{A}_{16}$ as $\tau \to 0$ and to
match $J_1$ as $\tau \to \infty$, using a $\tau$-dependent
modulation of the correction term.

The precise algebraic form of $\mathcal{A}_{16}$ and
$\mathcal{A}_{17}$ printed in Evans (1998) Eqs. (A16) and (A17)
does not, in the author's reading, reduce to the Taylor expansion
{eq}`eq-I-taylor` in the constant-extinction limit at second
order — a small inconsistency that may be a typographical issue
in the published formulas. We do not reproduce these formulas
here; see Evans (1998) pages 443–444 for the printed expressions.
For the purposes of this review, {eq}`eq-A15`–{eq}`eq-I-taylor`
establish the *exact* cell-integration behaviour in the two
limits (constant-$k$; small-$\tau$), and any efficient
approximation used by the SHDOM implementation can be verified
against these.

(app-shdom-cellint-face-interp)=
## Entering-radiance interpolation

The integration described above assumes that the entering
radiance $I(0)$ at the ray's intersection with the upstream cell
face is known. In general, the ray enters the cell face at a
point that is not a grid node, so $I(0)$ is obtained by bilinear
interpolation of the four corner-grid-point radiance values that
bracket the entry point on the face. Evans (1998) notes that this
interpolation is systematic (deterministic rather than stochastic
as in Monte Carlo methods) and, together with the source-function
integration error within the cell, is the dominant contributor to
SHDOM's spatial truncation error. The adaptive-grid splitting
criterion derived in [](adaptive-grid.md) targets the
source-function component of this error; residual interpolation
error, by contrast, is controlled only by the *base* grid
resolution (Evans, 1998, §2), which is why the SHDOM user must
supply a base grid fine enough to resolve the radiance field, not
only the source function.
