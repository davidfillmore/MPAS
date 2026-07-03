(sec-terrain)=
# Terrain-following meshes via altitude-grid remap

The flat-terrain assumption of [](discretization.md)
($z_{\mathrm{g}}(\boldsymbol{x}_h) = 0$) is lifted not by adding
cross-derivative metric terms to the operator, but by solving on a
uniform altitude grid and remapping to and from the MPAS
terrain-following coordinate. This keeps the assembled operator
identical in form to {eq}`eq-L-full` — symmetric, two points per face,
SPD after negation — and reuses the entire PCG solver unchanged.

## Altitude grid and cut cells

A single uniform vertical grid is built per domain with $K$ bands
bounded by the faces
$z^{\mathrm{f}}_k = z_{\mathrm{b}} + (k-1)\,\Delta z_{\mathrm{P}}$,
$k = 1, \dots, K{+}1$, with
$\Delta z_{\mathrm{P}} = (z_{\mathrm{t}} - z_{\mathrm{b}})/K$, where the
extents $z_{\mathrm{b}} = \min_i z_{\mathrm{g}}(\boldsymbol{x}_i)$ and
$z_{\mathrm{t}}$ are reduced globally across MPI ranks so every rank
shares the same grid. Terrain detection is likewise a global reduction,
so a rank whose partition is locally flat still follows the terrain path
and the collective solve cannot deadlock. Each column $i$ is clipped by
its ground height $z_{\mathrm{g}}(\boldsymbol{x}_i)$: a cell is active
where its clipped air band has positive thickness, and the first active
band reaches down to the terrain,
$[\,z_{\mathrm{g}}(\boldsymbol{x}_i),\,z^{\mathrm{f}}_{k_1+1}]$, so that
the air column is conserved exactly,
$\sum_k \Delta z^{\mathrm{air}}_{k,i} = z_{\mathrm{t}} - z_{\mathrm{g}}(\boldsymbol{x}_i)$.
A thin first band (below $0.1\,\Delta z_{\mathrm{P}}$) is merged upward
to bound the cut-cell aspect ratio. The unknown of a cut cell sits at
its air centroid, so the ground Dirichlet distance is always half the
air thickness — the classical cut-cell half-thickness, reducing to
$\Delta z_{\mathrm{P}}/2$ on flat columns.

## Shaved faces and grounded walls

Lateral faces are shaved by the two adjacent ground heights. The
air–air overlap
$t_{\mathrm{open}} = \max(0,\,z^{\mathrm{f}}_{k+1} - \max(z^{\mathrm{f}}_k, z_{\mathrm{g},i}, z_{\mathrm{g},j}))$
carries the usual flux weight
$\varepsilon\,(\ell_e/d_e)\,t_{\mathrm{open}}$, which is face-intrinsic
and therefore symmetric. Where air abuts rock across the edge, the
overlap is a grounded Dirichlet wall at distance $d_e/2$, contributing
$\varepsilon\,\ell_e\,t_{\mathrm{wall}}/(d_e/2)$ to the cell's diagonal
and, paired with it, the same datum $\varphi = \varphi_{\mathrm{g}}$ to
the right-hand side. Cells below the terrain are inactive identity rows
($\varphi = 0$), a grounded interior, and the bottom face of the first
active cell is Dirichlet at the exact terrain height. Every wall and
bottom face thus carries a single equipotential datum: on a grounded
conductor with $\rho = 0$ the discrete constant state
$\varphi \equiv \varphi_{\mathrm{g}}$ satisfies
$A\boldsymbol{x} = \boldsymbol{b}$ exactly. The operator remains SPD, and
flat meshes — run through the identical builders with the native
(possibly stretched) interface profile — reproduce {eq}`eq-L-full` to
round-off, so a single geometry path serves both regimes.

## Conservative charge remap and field return

The charge density is transferred from the MPAS terrain-following layers
to the altitude bands by piecewise-constant overlap integration, which
preserves the column-integrated charge to machine precision. After the
solve, $\varphi$ and the cell-centered $\boldsymbol{E}$ are interpolated
back to the MPAS layer midpoints for output; the edge-normal component
$E_{\mathrm{n}}$ is retained on the altitude grid, whose per-column
heights are recorded by the `poisson_zmid` diagnostic.

## Manufactured-solution verification

A terrain manufactured solution over a cosine hill
$s(x, y) = h_0 \cos(k_x x) \cos(k_y y)$, with
$u = (z - s)/(z_{\mathrm{t}} - s)$ and
$\varphi = \sin(\tfrac{\pi}{2} u)$ (grounded on the terrain,
zero-Neumann at the top), is solved by the full Fortran pipeline under
combined horizontal/vertical refinement
$(15\,\mathrm{km}, K{=}25) \to (7.5\,\mathrm{km}, K{=}50) \to (3.75\,\mathrm{km}, K{=}100)$.
For a gentle $1\,\mathrm{km}$ hill the interior relative-$L^2$ error
converges at finest-two slope $1.978$ (second order, matching the flat
Cartesian rate) with final PCG residuals below $10^{-8}$. A steep
$6\,\mathrm{km}$ hill converges at $0.921$: the near-surface field over
strong slopes is first order, because the cut-cell reconstruction does
not enforce the vanishing of the discrete tangential electric field
along the sloped equipotential wall. This is the one residual
approximation of the scheme, and it is largest where the terrain is
steepest.

## Coupled terrain supercell

The supercell example of [](end-to-end-supercell.md) is repeated over a
smooth $1\,\mathrm{km}$ cosine hill of $20\,\mathrm{km}$ half-width.
Imposing terrain requires the host's `zgrid`-derived vertical metrics to
be recomputed consistently with the hill, so the dynamical core evolves
on a mesh that actually carries the terrain rather than on flat-plane
metrics. With the same $\alpha = \beta = 10^{-9}$ charging as the flat
example, the diagnostic sequence
({numref}`fig-terr-030`–{numref}`fig-terr-120`, a $10\,\mathrm{km}$
half-width storm-core mean) reaches a maximum potential of
$207\,\mathrm{MV}$ and peak field $4.5 \times 10^4\,\mathrm{V\,m^{-1}}$
at $90\,\mathrm{min}$, with charge density in
$[-0.13,\,0.65]\,\mathrm{nC\,m^{-3}}$ and final PCG residuals below
$10^{-10}$ throughout — magnitudes consistent with the flat supercell
run, as expected for one-way linear charging. The solid ground conceals
the sub-surface region, where the operator grounds the potential; the
storm develops above the hill and the field structure closes on the
grounded terrain, confirming that the conservative remap, cut-cell
operator, and field return compose correctly with the terrain-following
host dynamics.

```{figure} figures/tier_E_terrain_charge_coupled_030min.png
:name: fig-terr-030
:alt: Terrain-following charge-coupled supercell at 30 min.
:width: 100%

Terrain-following charge-coupled supercell over a $1\,\mathrm{km}$
cosine hill at 30 min, shown as a $10\,\mathrm{km}$ half-width mean
through the storm core. The left panel shows liquid water content with
negative (red) and positive (black/gray) charge-density contours; the
right panel shows electrostatic potential with electric-field
streamlines. The solid fill is the grounded terrain, below which the
operator holds $\varphi = 0$.
```

```{figure} figures/tier_E_terrain_charge_coupled_060min.png
:name: fig-terr-060
:alt: Terrain-following charge-coupled supercell at 60 min.
:width: 100%

As in {numref}`fig-terr-030`, but at 60 min, with maximum potential
$113\,\mathrm{MV}$ and peak field $2.6 \times 10^4\,\mathrm{V\,m^{-1}}$.
```

```{figure} figures/tier_E_terrain_charge_coupled_090min.png
:name: fig-terr-090
:alt: Terrain-following charge-coupled supercell at 90 min.
:width: 100%

As in {numref}`fig-terr-030`, but at 90 min. This frame has the largest
potential ($207\,\mathrm{MV}$) and field
($4.5 \times 10^4\,\mathrm{V\,m^{-1}}$) of the plotted terrain sequence.
```

```{figure} figures/tier_E_terrain_charge_coupled_120min.png
:name: fig-terr-120
:alt: Terrain-following charge-coupled supercell at 120 min.
:width: 100%

As in {numref}`fig-terr-030`, but at 120 min, after the core has cycled:
maximum potential $76\,\mathrm{MV}$ and field
$2.0 \times 10^4\,\mathrm{V\,m^{-1}}$. The remap, Poisson solve, and
field reconstruction remain coupled to the evolving convective state
over terrain.
```
