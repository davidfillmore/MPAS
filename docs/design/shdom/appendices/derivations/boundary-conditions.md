(app-shdom-bc)=
# Radiation boundary conditions

Radiative transfer in a finite-volume domain requires radiance
boundary conditions on all faces. SHDOM applies directional
inflow conditions on the top of the atmosphere, a
reflection-plus-emission bidirectional-reflectance distribution
function (BRDF) condition on the bottom, and either periodic wrap
or "independent scan" open conditions on the horizontal sides.
This section derives Eq. (A18) of Evans (1998) for the lower
boundary and specifies the top and horizontal-boundary
treatments.

(app-shdom-bc-top)=
## Top-of-atmosphere boundary: inflow radiance

At the top of the atmosphere the domain is bounded by empty
space; no diffuse radiation enters from above, but the collimated
solar beam does. In SHDOM's direct–diffuse decomposition
(see {ref}`app-shdom-sources-solar-motivation`), the
solar beam is integrated analytically along the solar ray by
Beer's law and appears in the SHDOM solver only as the
pseudo-source $S^{\odot}_{lm}$ of Eq. (A7). The *diffuse*
radiance boundary condition at the top is therefore usually

$$
  I^{\mathrm{diff}}(\mu, \phi) \;=\; 0
  \qquad \text{for all } \mu < 0 \text{ at } z = z_{\mathrm{top}},
$$ (eq-top-bc)

i.e. zero incoming diffuse radiation from above. An isotropic
diffuse-incident option is also supported (e.g. for incoming
cosmic-microwave-background-calibrated thermal radiation in
far-IR applications) by setting {eq}`eq-top-bc` to a non-zero
constant, but is not relevant to the present coupling work.

(app-shdom-bc-A18)=
## Lower boundary: BRDF reflection plus thermal emission: Eq. (A18)

At the lower surface ($z = z_{\mathrm{g}}$) the upward-going
radiance in direction $(\mu, \phi)$ with $\mu > 0$ results from
two physical mechanisms:

- Reflection of incident downward radiance $I(-\mu', \phi')$
  (with $\mu' \in (0, 1]$) by a surface with bidirectional-
  reflectance distribution function
  $\rho(\mu, \phi; -\mu', \phi')$. The reflected radiance into
  direction $(\mu, \phi)$ integrates over all incident downward
  directions weighted by $\mu'$ (projected solid angle).
- Thermal emission at the surface temperature $T_{\mathrm{s}}$
  with directional emissivity
  $\varepsilon(\mu, \phi) = 1 - \int \rho(\mu, \phi; -\mu', \phi')\,\mu'\,d\mu'\,d\phi'/\pi$
  (Kirchhoff's law); the thermal radiance is
  $\varepsilon(\mu, \phi)\,B(T_{\mathrm{s}})$.

Combining the two contributions, the upward boundary radiance is

$$
  \boxed{\;
    I(\mu, \phi)
    \;=\; \frac{1}{\pi}\!
    \int_{0}^{2\pi}\!\!\!\int_{0}^{1}
    \Bigl\{\rho(\mu, \phi; -\mu', \phi')\,I(-\mu', \phi')
    + [1 - \rho(\mu, \phi; -\mu', \phi')]\,B(T_{\mathrm{s}})\Bigr\}\,
    \mu'\,d\mu'\,d\phi',
  \;}
$$ (eq-A18)

which is Eq. (A18) of Evans (1998). The $(1/\pi)$ prefactor
absorbs the Lambertian normalization: for a Lambertian surface of
albedo $\alpha$ the BRDF is $\rho_{L} = \alpha/\pi$ (so
$\int \rho_{L} \mu'\,d\mu'\,d\phi' / \pi = \alpha$, consistent
with energy conservation).

(app-shdom-bc-special)=
## Lambertian and Fresnel specializations

**Lambertian surface.** Setting $\rho = \alpha/\pi$ in
{eq}`eq-A18` and using
$\int_{0}^{2\pi}\int_{0}^{1} \mu'\,d\mu'\,d\phi' = \pi$,

$$
  I(\mu, \phi) \;=\; \frac{\alpha}{\pi}\,F^{\downarrow}
                   \;+\; (1 - \alpha)\,B(T_{\mathrm{s}}),
$$ (eq-lambertian)

with $F^{\downarrow} \equiv \int_{0}^{2\pi}\int_{0}^{1}
I(-\mu', \phi')\,\mu'\,d\mu'\,d\phi'$ the downward irradiance on
the horizontal surface. This is the textbook Lambertian form:
reflected radiance is isotropic at rate
$(\alpha/\pi)\,F^{\downarrow}$ and thermal emission contributes
$(1 - \alpha) B(T_{\mathrm{s}})$, the "gray-body" emissivity
being $1 - \alpha$ by Kirchhoff's law for a Lambertian surface.

**Fresnel surface.** For a specular surface the BRDF collapses to
a $\delta$-function in angles,

$$
  \rho^{\mathrm{Fresnel}}(\mu, \phi; -\mu', \phi')
  \;=\; r(\mu)\,\frac{\delta(\mu - \mu')\,\delta(\phi - \phi' - \pi)}{\mu},
$$

with $r(\mu)$ the Fresnel reflectance for unpolarized light at
incidence angle $\arccos\mu$. Substituting into {eq}`eq-A18` and
integrating out the $\delta$-functions gives
$I(\mu, \phi) = r(\mu)\,I(-\mu, \phi + \pi) + [1 - r(\mu)]\,B(T_{\mathrm{s}})$,
the specular reflection plus thermal emission. Because the
specular reflection requires evaluating $I$ at a direction that
is not in general a discrete ordinate, SHDOM handles Fresnel
surfaces in separate routines rather than via the general
{eq}`eq-A18` path (Evans, 1998, appendix §d).

**General BRDF.** For an arbitrary tabulated BRDF
(Rahman–Pinty–Verstraete model, MODIS BRDF retrievals, etc.) the
integral in {eq}`eq-A18` is evaluated on the same
discrete-ordinate grid used for the volumetric solution; this is
consistent because the surface-reflected radiance lives on
exactly those ordinates.

(app-shdom-bc-horizontal)=
## Horizontal boundary conditions: periodic and open

In the horizontal, SHDOM supports two boundary treatments,
independently selectable in each of $X$ and $Y$:

- **Periodic wrap.** Ray segments exiting through one side
  re-enter at the opposite side; neighbour pointers in the grid
  data structure are set to wrap. Grid values on the two wrap
  faces are identical by construction. This is appropriate for
  idealized tests (e.g. doubly-periodic LES or supercell test
  cases).
- **Open ("independent scan").** Boundary columns are made
  *independent* of each other and of the interior, so each
  boundary column sees plane-parallel inflow from the sides with
  no horizontal coupling. Operationally the boundary cells point
  to themselves as neighbours rather than to interior cells, so
  the ray-tracing "sweep" through the grid treats the boundary
  column as if it were surrounded by identical copies of itself.
  For 2D ($X$–$Z$) transfer this reproduces plane-parallel
  radiance on the boundary columns.

Mixed (periodic in one direction, open in the other) is
permitted. The coupling narrative's supercell test case uses
periodic $\times$ periodic in $(X, Y)$; limited-area
open-boundary configurations are deferred to later phases of the
work.
