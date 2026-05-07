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

The lateral-edge values {eq}`eq-G-edge-E` determine the horizontal
or locally tangential part of the electric field. To reconstruct
that part at cell centres, we use MPAS-A's existing
radial-basis-function (RBF) machinery in
`mpas_rbf_interpolation.F`. This machinery takes a set of lateral
edge-normal values and, for each cell $i$, fits a vector
$\bE^{h}_i$ such that for each horizontal edge $e \in \partial i$:

$$
  \bE^{h}_i \cdot \bn_e \;\approx\; E_n^{(e, k)}.
$$ (eq-G-rbf)

The fit is a linear map built from the cell's edge geometry plus
a set of RBF coefficients precomputed at mesh init. The
coefficients are mesh-only and do not depend on $\varphi$: they
are computed once and stored. The reconstruction call at each
time step is therefore a matrix-vector product of fixed dimension
(12–24 for a typical hex) per cell. Because all source normals are
lateral-edge normals, this RBF reconstruction does not determine
the vertical component of $\bE$.

The vertical component is computed separately from the same
layer-centred potential used in the vertical operator:

$$
  E_{\mathrm{vert}}^{(i,k)}
  \;=\; -\partial_z \varphi
  \;\approx\;
  -\frac{\varphi_{i,k+1} - \varphi_{i,k-1}}
         {z_{i,k+1} - z_{i,k-1}}
  \quad \text{for interior } k,
$$ (eq-G-vertical-E)

with one-sided, boundary-condition-consistent differences at the
ground and model top. At the ground face the prescribed
Dirichlet value is used; at the model top the Neumann condition
sets the upper-face vertical field to zero.

The final cell-centred vector is therefore

$$
  \bE_i = \bE^{h}_i
          + E_{\mathrm{vert}}^{(i,k)} \, \hat{\boldsymbol{z}}_i,
$$ (eq-G-full-E)

where $\hat{\boldsymbol{z}}_i$ is the Cartesian vertical unit
vector on a flat periodic mesh and the local radial unit vector on
a spherical mesh.

The same machinery is used to reconstruct cell-centred velocity
vectors from edge-normal momentum in the MPAS-A dynamics core;
full details are given in
{cite:t}`skamarock2012multiscale,ringler2010unified`.

## Output to Registry

The reconstructed vector is written to the Registry field
`E_vector` with dimensions `(R3, nVertLevels, nCells)`. The
vector combines the RBF-reconstructed horizontal component with
the finite-difference vertical component and is expressed in the
model's Cartesian coordinate system
($x, y, z$ ECEF-style on the sphere, or flat Cartesian for
doubly-periodic configurations). A user may derive the magnitude
$|\bE|$ from `E_vector` in post-processing.
