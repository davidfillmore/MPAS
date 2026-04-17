(app-shdom-rte)=
# Radiative transfer equation: differential and integral forms

This section derives both the differential form of the monochromatic
radiative transfer equation (Eq. (A11) of Evans, 1998) from a
photon-number balance, and its integration along a ray to yield the
integral form (Eq. (2) of Evans, 1998) used by SHDOM to advance the
radiance. It ends by placing the source-function iteration in its
operator-theoretic setting as a Picard ("$\Lambda$") iteration.

(app-shdom-rte-diff)=
## Differential form from photon balance

Consider a pencil of monochromatic rays of infinitesimal cross-section
traveling in a unit direction $\hat{\Omega}$ through the atmosphere.
Let $I(s, \hat{\Omega})$ denote the radiance at position $s$ along the
ray, where $s$ is arclength measured from some reference point on the
ray (we will take $s = 0$ at a boundary or upstream grid point). Over
an infinitesimal step $ds$, photons are removed from the beam by
extinction at the rate

$$
  -k(s)\,I(s, \hat{\Omega})\,ds,
$$ (eq-rte-extinction-loss)

where $k(s)$ is the volume extinction coefficient (absorption plus
scattering, with dimension $[L]^{-1}$). Photons are added by
in-scattering from all other directions and by thermal emission.
Writing the combined source per unit path length as
$k(s)\,J(s, \hat{\Omega})$, where $J$ is the *source function* with
the same units as radiance, the photon-number balance is

$$
  dI(s, \hat{\Omega})
  = -k(s)\,I(s, \hat{\Omega})\,ds + k(s)\,J(s, \hat{\Omega})\,ds.
$$ (eq-rte-balance)

Dividing by $ds$ gives the differential form of the radiative transfer
equation,

$$
  \boxed{\;
    \frac{dI(s, \hat{\Omega})}{ds}
    = -k(s)\,\bigl[I(s, \hat{\Omega}) - J(s, \hat{\Omega})\bigr],
  \;}
$$ (eq-A11)

reproducing Eq. (A11) of Evans (1998).

The source function comprises three physical processes:

$$
  J(s, \hat{\Omega})
  \;=\; \underbrace{\frac{\omega(s)}{4\pi}
      \int_{4\pi} P(s; \hat{\Omega}, \hat{\Omega}')\,
      I(s, \hat{\Omega}')\,d\hat{\Omega}'}_{\text{in-scattering}}
   \;+\; \underbrace{(1-\omega(s))\,B\!\bigl(T(s)\bigr)}_{\text{thermal emission}}
   \;+\; \underbrace{S^{\odot}(s, \hat{\Omega})}_{\text{solar pseudo-source}},
$$ (eq-J-physical)

in which $\omega(s)$ is the single-scattering albedo,
$P(\hat{\Omega}, \hat{\Omega}')$ is the scattering phase function
normalized so that $\int_{4\pi} P\,d\hat{\Omega}/(4\pi) = 1$, $B$ is
the Planck function at temperature $T(s)$, and $S^{\odot}$ is the
singly-scattered direct-beam pseudo-source. Legendre expansions of
the scattering integral are carried out in [](scattering.md); the
thermal and solar pieces are derived in [](sources.md).

(app-shdom-rte-int)=
## Integral form via the integrating factor

Equation {eq}`eq-A11` is a first-order linear inhomogeneous ODE along
the ray. Rearranging,

$$
  \frac{dI}{ds} + k(s)\,I(s) = k(s)\,J(s),
$$ (eq-rte-linODE)

in which we suppress the direction argument $\hat{\Omega}$ for
brevity. Define the optical depth between two points $a$ and $b$ on
the ray,

$$
  \tau(a, b) \equiv \int_{a}^{b} k(s')\,ds',
$$ (eq-tau-def)

and the integrating factor $\mu(s) \equiv e^{\tau(0, s)}$ (not to be
confused with the zenith cosine $\mu$ introduced later). Multiplying
{eq}`eq-rte-linODE` by $\mu(s)$ yields

$$
  \frac{d}{ds}\!\left[e^{\tau(0, s)}\,I(s)\right]
  = e^{\tau(0, s)}\,k(s)\,J(s),
$$ (eq-rte-integrable)

which we integrate from $s' = 0$ to $s' = s$:

$$
  e^{\tau(0, s)}\,I(s) - I(0)
  = \int_{0}^{s} e^{\tau(0, s')}\,k(s')\,J(s')\,ds'.
$$ (eq-rte-integrated-raw)

Multiplying both sides by $e^{-\tau(0, s)}$ and using the identity
$\tau(0, s) - \tau(0, s') = \tau(s', s)$ (which follows directly from
{eq}`eq-tau-def` for $0 \le s' \le s$):

$$
  \boxed{\;
    I(s) \;=\; e^{-\tau(0, s)}\,I(0)
    \;+\; \int_{0}^{s} e^{-\tau(s', s)}\,k(s')\,J(s')\,ds',
  \;}
$$ (eq-eq2)

which reproduces Eq. (2) of Evans (1998). The first term transmits
the boundary-entering radiance $I(0)$ through the entire path optical
depth $\tau(0, s)$; the second term accumulates source contributions
from every point $s' \in [0, s]$ along the ray, each attenuated by
the *remaining* path optical depth $\tau(s', s)$.

Two remarks on Evans's notation:

- Evans writes the first-term exponential as
  $\exp\!\bigl[-\!\int_{0}^{s}\!k(s'')\,ds''\bigr]$, i.e. our
  $e^{-\tau(0, s)}$; and the second-term exponential as
  $\exp\!\bigl[-\!\int_{s'}^{s}\!k(t)\,dt\bigr]$, i.e. our
  $e^{-\tau(s', s)}$. The two forms are identical; we use the
  $\tau(\cdot, \cdot)$ notation for compactness in later appendices.
- The appearance of $k(s')$ inside the $J$-integral is *essential*:
  it reflects that $J$ is the emission/scattering source *per unit
  extinction event*, not per unit path length. Equivalently, the
  product $(Jk)$ is the emission coefficient per unit path length,
  which is the physically natural quantity; accordingly SHDOM treats
  $(Jk)$ — not $J$ alone — as linear in $s$ during cell integration
  (see [](cell-integration.md)).

(app-shdom-rte-picard)=
## Picard iteration (the "$\Lambda$" operator)

For fixed inflow boundary-entering radiance $I(0)$,
equation {eq}`eq-eq2` is a linear integral mapping from the source
function $J(\cdot, \hat{\Omega})$ along the ray to the outgoing
radiance $I(s, \hat{\Omega})$. Sweeping every discrete ordinate
through every cell of the computational grid composes these
per-ray mappings into a global linear operator $\mathcal{L}$ acting
on the source-function field $J(\mathbf{r}, \hat{\Omega})$:

$$
  I \;=\; \mathcal{L}\,J \;+\; I_{\mathrm{bc}},
$$ (eq-rte-Lop)

where $I_{\mathrm{bc}}$ collects the contribution of the transmitted
boundary inflow (the first term in {eq}`eq-eq2`, summed over all
rays entering the domain through the boundaries).

Substituting {eq}`eq-rte-Lop` into the scattering contribution to
{eq}`eq-J-physical` rearranges the problem into a linear equation
for $J$ alone,

$$
  J \;=\; \mathcal{S}\,I + Q
  \;=\; \mathcal{S}\mathcal{L}\,J + \bigl(\mathcal{S}\,I_{\mathrm{bc}} + Q\bigr)
  \;\equiv\; \mathcal{K}\,J + f,
$$ (eq-rte-Kop)

where $\mathcal{S}$ is the angular integral against the phase
function (evaluated efficiently in spherical-harmonic space, see
[](scattering.md)), $Q$ collects the thermal and direct-solar pieces
of {eq}`eq-J-physical`, and $\mathcal{K} \equiv \mathcal{S}\mathcal{L}$
is the composite "scatter-then-stream" operator. In the
astrophysical literature $\mathcal{L}$ is known as the $\Lambda$
operator and {eq}`eq-rte-Kop` as the $\Lambda$ (or Picard) equation
{cite:p}`stenholm1991lambda`.

SHDOM solves {eq}`eq-rte-Kop` by Picard iteration,

$$
  J^{(n+1)} \;=\; \mathcal{K}\,J^{(n)} + f,
$$ (eq-rte-Picard)

which, in implementation, is the four-step cycle sketched in §2 of
the main text of Evans (1998) and recapitulated here:

1. *Transform* the spherical-harmonic source $J_{lm}^{(n)}$ to
   discrete ordinates $J_{jk}^{(n)}$ via Eq. (A2) of Evans (1998);
   this is the forward transform of
   [](spherical-harmonics.md).
2. *Stream*: integrate {eq}`eq-eq2` along every discrete ordinate
   through every cell to obtain $I_{jk}^{(n+1)}$; the
   cell-integration formula is derived in
   [](cell-integration.md).
3. *Transform* $I_{jk}^{(n+1)}$ back to spherical harmonics
   $I_{lm}^{(n+1)}$ via Eq. (A3) of Evans (1998); this is the
   inverse transform of [](spherical-harmonics.md).
4. *Evaluate scattering*: compute $J_{lm}^{(n+1)}$ from
   $I_{lm}^{(n+1)}$ via Eq. (A5) of Evans (1998), which in
   spherical harmonic space is pointwise multiplication by
   $\omega\chi_{l}/(2l+1)$, derived in
   [](scattering.md).

Iteration is terminated when the RMS change in $J$ between
successive iterates falls below a prescribed solution criterion. The
spectral radius of $\mathcal{K}$ approaches unity for conservative
scattering ($\omega \to 1$) in optically thick media, causing the
iteration to stall; the sequence-acceleration scheme
([](acceleration.md)) addresses this by geometric-series
extrapolation of the Picard iterates.
