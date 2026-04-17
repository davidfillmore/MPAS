(app-E-recon)=
# Edge-normal and cell-centred E-field reconstruction

## Edge-normal component

From $\bE = -\nabla \varphi$ and the edge-normal finite-difference
gradient {eq}`eq-A-face-grad`, the edge-normal $\bE$ component
at edge $e$ between cells $i$ and $j = \mathrm{nbr}(i, e)$ is

$$
  E_n^{(e, k)} \;=\;
  -(\nabla \varphi) \cdot \bn_e
  \;\approx\;
  -\frac{\varphi_{j, k} - \varphi_{i, k}}{\dcen}.
$$ (eq-G-edge-E)

This is stored in the Registry field `E_normal`.

## Cell-centred vector reconstruction

To reconstruct the full vector $\bE$ at cell centres, we use
MPAS-A's existing radial-basis-function (RBF) machinery in
`mpas_rbf_interpolation.F`. This machinery takes a set of
edge-normal values and, for each cell $i$, fits a 3-component
vector $\bE_i$ such that for each edge $e \in \partial i$:

$$
  \bE_i \cdot \bn_e \;\approx\; E_n^{(e, k)}.
$$ (eq-G-rbf)

The fit is a linear map built from the cell's edge geometry plus
a set of RBF coefficients precomputed at mesh init. The
coefficients are mesh-only and do not depend on $\varphi$: they
are computed once and stored. The reconstruction call at each
time step is therefore a matrix-vector product of fixed dimension
(12–24 for a typical hex) per cell.

The same machinery is used to reconstruct cell-centred velocity
vectors from edge-normal momentum in the MPAS-A dynamics core;
full details are given in
{cite:t}`skamarock2012multiscale,ringler2010unified`.

## Output to Registry

The reconstructed vector is written to the Registry field
`E_vector` with dimensions `(R3, nVertLevels, nCells)`. The
vector is expressed in the model's Cartesian coordinate system
($x, y, z$ ECEF-style on the sphere, or flat Cartesian for
doubly-periodic configurations). A user may derive the magnitude
$|\bE|$ from `E_vector` in post-processing.
