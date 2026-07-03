(sec-conclusions)=
# Conclusions

This narrative has described the first elliptic Poisson solver for the
MPAS-A unstructured Voronoi mesh, solving
$\nabla \cdot (\varepsilon \nabla \varphi) = -\rho$ for electrostatic
applications. The discrete operator reuses MPAS-A's existing C-grid mesh
metrics, is symmetric with a positive-definite negation after volume
integration, and is solved by preconditioned conjugate gradients with no
external linear-algebra dependency. The results reported here establish
second-order horizontal convergence in the Cartesian MMS mode and
document a localized pentagon-defect limitation for the spherical MMS. A
controlled mesh-sensitivity study shows that the second-order accuracy
is confined to near-uniform hexagonal meshes, collapsing to a
near-stalled $\sim\!0.35$ rate on the irregular, general
centroidal-Voronoi tessellations of variable-resolution configurations,
without any degenerate cells. Idealized point-charge and tripole
diagnostics demonstrate physically interpretable fields, and the one-way
charge-coupled supercell diagnostic path produces time-sampled $\rho$,
$\varphi$, and $\boldsymbol{E}$ fields with tightly converged PCG
residuals, including over terrain-following meshes handled by a
conservative altitude-grid remap ([](terrain-following.md)) that
preserves the SPD operator and the second-order interior rate.

The principal limitations are explicit: the near-surface field over
steep terrain is first order; restoring second order on irregular and
variable-resolution meshes requires a higher-order or finite-element
reconstruction; and large-scale parallel scaling, physical
electrification, lightning discharge, and two-way feedback are outside
the present scope. The capability thus established is the foundation for
future work on lightning parameterizations, lightning NO$_x$ emission,
and a global-atmospheric-electric-circuit solver that extends the same
operator with variable conductivity.
