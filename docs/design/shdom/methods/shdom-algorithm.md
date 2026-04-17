(sec-shdom-algorithm)=
# SHDOM as a 3D RT engine

The Spherical Harmonic Discrete Ordinate Method
{cite:p}`evans1998shdom` solves the monochromatic radiative transfer
equation in three dimensions on a Cartesian grid by carrying two
complementary representations of the radiance field and converting
between them at each step of a Picard iteration. The radiance
itself is advanced along discrete ordinates by integration of the
integral form of the transfer equation across grid cells; the
source function, which couples the ordinates through in-scattering,
is represented in spherical harmonics on each grid point.
Transforming between the two representations through a
reduced-Gaussian discrete-ordinate grid permits the scattering
integral, which is a local multiplication in spherical harmonic
space, to be evaluated in $O(N^{3/2})$ operations per grid point
rather than the $O(N^2)$ that a direct discrete-ordinate
convolution would require, where $N$ is the number of discrete
ordinates. Full step-by-step derivations of the spherical-harmonic
representation, the scattering integral, the cell-integration
algorithm, and the source-function iteration follow the notation
of {cite:t}`evans1998shdom`; they will appear as appendix pages of
this narrative in a later revision.

SHDOM maintains an adaptive Cartesian cell structure that
concentrates degrees of freedom on regions of strong radiative
gradient; resolves the forward-peaked phase function by delta-M
truncation {cite:p}`wiscombe1977deltam` and the TMS single-scatter
correction of {cite:t}`nakajima1988tms`; supports periodic or
collimated-inflow upper boundary conditions and arbitrary
bidirectional surface reflectances at the lower boundary; and
accelerates the Picard ("$\Lambda$") source-function iteration by
geometric-sequence extrapolation {cite:p}`stenholm1991lambda`.
Radiometric outputs include hemispheric fluxes, mean radiance, net
flux convergence, and radiance at specified viewing directions.

Three features of the SHDOM implementation as it is embedded here
bear emphasizing, because they delimit the scope of this narrative
rather than of SHDOM as a method. First, the geometric domain is
Cartesian and rectilinear — SHDOM does not know about MPAS-A's
Voronoi mesh and does not need to; the interface between the two is
the mesh mapper of [](mesh-mapping.md). Second, the phase function
in this Phase 1 configuration is reduced to its Henyey–Greenstein
form {cite:p}`henyey1941`, parameterized by asymmetry $g$ alone;
extension to fully tabulated Legendre-series phase functions is
supported by SHDOM but unused here. Third, the spectral dimension
is pre-compiled at the band level: Phase 1 operates in a single
visible band, and extension to multiple bands is a scoped future
phase (see [](../discussion.md)). Polarization, which SHDOM
supports as an option, is unused throughout.
