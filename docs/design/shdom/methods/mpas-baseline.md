(sec-mpas-baseline)=
# MPAS-A baseline radiation pipeline

MPAS-A solves the nonhydrostatic equations on a centroidal Voronoi
tessellation of the sphere with C-grid staggering of the prognostic
variables {cite:p}`skamarock2012mpas`. Physical parameterizations
are dispatched from a per-timestep driver that traverses the cells
of each MPI rank's local sub-mesh and calls the selected
microphysics, boundary-layer, convection, and radiation schemes
column by column. For shortwave and longwave radiation, the bundled
choice in the present release is the Rapid Radiative Transfer Model
for GCMs (RRTMG) of {cite:t}`mlawer1997rrtm` and
{cite:t}`iacono2008rrtmg`, which integrates the radiative transfer
equation along each column independently using a correlated-$k$
discretization of the absorption spectrum and a Monte Carlo
independent-column approximation for sub-grid cloud overlap
{cite:p}`pincus2003mcica,barker2003monte`. Cloud optical properties
are reconstructed per $g$-point on a stochastic sample of
binary-cloud sub-columns; the resulting per-column heating rates
are tendered back to the dynamical core through the standard
physics tendency fields.

The pipeline is column-wise by construction: each MPAS-A cell's
radiation solve is independent of every other cell, there is no
horizontal photon transport, and the only coupling between columns
is mediated indirectly by the advected cloud field. The hybrid
1D/3D framework described in this narrative preserves this baseline
unchanged and inserts an additional, optional 3D pass alongside it;
when the framework is disabled the radiation pipeline reduces
exactly to the existing RRTMG call sequence.
