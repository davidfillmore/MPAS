(tier-d)=
# Tier D — Parallel performance

Tier D reports the solver's parallel scaling behaviour under the
two non-trivial preconditioners of [](../solver.md). Two
complementary measurements are planned.

## Strong scaling

A single global $30$ km quasi-uniform MPAS-A mesh ($\sim 65{,}000$
cells $\times$ 50 vertical layers, for approximately
$3.3 \times 10^{6}$ unknowns) is fixed, and the same problem is
solved at rank counts $\mathrm{nproc} \in \{16, 64, 256, 1024\}$.
Reported quantities per run are: wall time per PCG solve; PCG
iteration count at convergence; and the MPI-communication
fraction of the per-iteration wall time, extracted from on-line
timers around the matvec halo exchange and the two inner-product
global reductions.

## Weak scaling

A sequence of meshes of decreasing nominal spacing
$h \in \{60, 30, 15\}$ km is paired with matched rank counts
$\{64, 256, 1024\}$ so that the ratio of cells per rank stays
approximately fixed at $\sim 10^3$. Reported quantities per run
are the same as for the strong-scaling study.

Both studies are repeated for the Jacobi and block symmetric
Gauss–Seidel preconditioners. The quantity of primary interest is
how the PCG iteration count grows with problem size at fixed
preconditioner: unpreconditioned or diagonally-preconditioned CG
on a 3D Poisson operator is expected to scale as
$\mathcal{O}(h^{-1})$ in iteration count, and the measured growth
rate will motivate — or not — the Phase 2 replacement of the
in-tree PCG by an algebraic multigrid method such as hypre's
BoomerAMG {cite:p}`falgout2006hypre`. No pass/fail criterion is
applied to Tier D; the reported scaling curves are the output.
