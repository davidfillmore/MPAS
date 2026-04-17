(app-shdom-adaptive)=
# Adaptive-grid cell-splitting criterion

SHDOM's adaptive grid adds spatial resolution where the source
function varies rapidly and the cell-integration formulas of
[](cell-integration.md) therefore incur the largest error. The
criterion for splitting a cell is a *local* measure of
source-function variability across the cell, evaluated
directionally-averaged on the sphere and modulated by the cell's
optical thickness. This section derives that criterion,
reproducing Eqs. (A8)–(A10) of Evans (1998).

(app-shdom-adaptive-motivation)=
## Source-function difference as error proxy

The cell-integration formulas {eq}`eq-A15` and {eq}`eq-I-taylor`
evaluate

$$
  \mathcal{I}_{\mathrm{cell}}
  \;\equiv\; \int_{0}^{s} e^{-\tau(s', s)}\,(Jk)(s')\,ds',
$$ (eq-Icell-def)

which approximates the contribution of the cell to the outgoing
radiance under the linear-$(Jk)$ parameterization {eq}`eq-A14`.
The truncation error of {eq}`eq-Icell-def` is controlled by the
*departure* of the true $(Jk)$ profile from linear, which at
leading order in cell size is proportional to the second
derivative $\partial^{2}(Jk)/\partial s^{2}$ (standard
trapezoidal-rule truncation). A finite-difference proxy for this
second derivative across a cell bounded by grid points $1$ and
$2$ is the first difference $k_2\,J^{(2)} - k_1\,J^{(1)}$ — the
change in $(Jk)$ across the cell — which is direction-dependent
because $J$ depends on $\hat{\Omega}$.

(app-shdom-adaptive-A8)=
## Angle-averaged difference: Eq. (A8)

Because a single cell-splitting decision must be made for a cell
carrying *all* discrete ordinates, the direction-dependent
difference is reduced to a single scalar by taking the RMS over
the sphere. Define

$$
  \boxed{\;
    |\Delta J|
    \;\equiv\; \frac{1}{\bar{k}}\,\Biggl\{\frac{1}{4\pi}
      \int_{0}^{2\pi}\!\!\!\int_{-1}^{1}
      \bigl[k_2\,J^{(2)}(\mu, \phi) - k_1\,J^{(1)}(\mu, \phi)\bigr]^{2}\,
      d\mu\,d\phi\Biggr\}^{1/2},
  \;}
$$ (eq-A8)

which is Eq. (A8) of Evans (1998). The factor $1/\bar{k}$ (with
$\bar{k}$ the cell-averaged extinction) converts the $(Jk)$
difference into units of source function, making $|\Delta J|$
dimensionally comparable to $J$ itself. The factor $1/(4\pi)$ is
the reciprocal of the unit-sphere solid angle; $|\Delta J|$ is
therefore a directional RMS average of the "scaled source-function
difference" across the cell.

(app-shdom-adaptive-A9)=
## Evaluation in spherical-harmonic space: Eq. (A9)

Because SHDOM already carries the source function in spherical
harmonic representation as part of the Picard iteration,
converting {eq}`eq-A8` to a sum over $(l, m)$ modes is
computationally advantageous. Applying Parseval's identity —
which for the orthonormal basis of
{ref}`app-shdom-sphharm-Ylm` reads

$$
  \int_{0}^{2\pi}\!\!\int_{-1}^{1} [f(\mu, \phi)]^{2}\,d\mu\,d\phi
  \;=\; \sum_{l, m} f_{lm}^{2},
$$ (eq-Parseval)

with $f_{lm}$ the spherical-harmonic coefficients
{eq}`eq-I-lm-continuous` — to the integrand of {eq}`eq-A8`,

$$
  \boxed{\;
    |\Delta J|
    \;=\; \frac{1}{\bar{k}}\,
          \Biggl\{\frac{1}{4\pi}
          \sum_{l, m}
          \bigl[k_2\,J_{lm}^{(2)} - k_1\,J_{lm}^{(1)}\bigr]^{2}
          \Biggr\}^{1/2},
  \;}
$$ (eq-A9)

reproducing Eq. (A9) of Evans (1998). This is the form used at
runtime: the sum has $O(N_{\mathrm{lm}}) = O(N)$ complexity versus
the $O(N)$ complexity of {eq}`eq-A8` on the discrete-ordinate
grid, but avoids converting $J$ back to angular space solely for
the splitting decision.

(app-shdom-adaptive-A10)=
## Optical-thickness modulation and the splitting criterion: Eq. (A10)

The directional source-function difference $|\Delta J|$ is a pure
*property-field* quantity; it does not yet account for how much
the cell's contribution weights into the radiance at its
boundary. Two limiting cases illustrate:

- **Optically thin** ($\tau \to 0$): the cell contributes
  $\mathcal{I}_{\mathrm{cell}} \to s\,(J_0 k_0 + J_1 k_1)/2 \to
  0$ regardless of $|\Delta J|$; the radiance at the outgoing
  face differs from the entering face only by a small amount, so
  the cell-integration error is small even if the source function
  changes sharply.
- **Optically thick** ($\tau \to \infty$):
  $\mathcal{I}_{\mathrm{cell}} \to J_1 / k$ (saturation), and the
  entering radiance is fully absorbed; again cell-integration
  error stays bounded even for large $|\Delta J|$, because
  $\tau \to \infty$ saturates the weight of any interior
  source-function variation.

Between these limits, the cell's radiance-level contribution
carries a weight proportional to $(1 - e^{-\tau})$ (the
"transmission deficit" factor that appeared in {eq}`eq-A15`),
which is zero at $\tau = 0$ and saturates at $1$ as
$\tau \to \infty$. The cell-splitting criterion is accordingly

$$
  \boxed{\;
    C \;\equiv\; |\Delta J|\,\bigl(1 - e^{-\tau}\bigr),
    \qquad \tau = \bar{k}\,d,
  \;}
$$ (eq-A10)

reproducing Eq. (A10) of Evans (1998), where $d$ is the cell
width along the splitting axis and $\bar{k}$ is the average cell
extinction. Cells with $C$ above a user-specified "splitting
accuracy" threshold are marked for subdivision. In 3D the
criterion is computed on each of the three candidate splitting
axes; the axis with the largest $C$ is chosen. The criterion is
averaged over the four edges of a cell crossing a potential
splitting plane (Evans, 1998, appendix §b).

(app-shdom-adaptive-tree)=
## Tree-structured adaptive refinement

The adaptive grid is organized as a binary tree rooted at each
base cell. Only the leaves (end nodes) participate in the
radiative transfer computation; interior tree nodes retain
parent–child pointers for cell-retrieval during ray tracing. At
each solution iteration, leaves with the highest $C$ values
(above threshold) are split; the cell-splitting threshold is
gradually lowered over iterations until the user-specified target
is reached, producing a progressively finer grid rather than a
one-shot refinement. Grid smoothness is maintained by splitting
additional cells in neighbourhoods of highly refined regions, to
ensure no two adjacent cells differ by more than one refinement
level {cite:p}`evans1998shdom`. Ray-tracing through the tree
uses neighbour pointers that are refreshed after each splitting
pass.
