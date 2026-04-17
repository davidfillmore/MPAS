(sec-mpi)=
# MPI redistribution between heterogeneous decompositions

MPAS-A partitions its Voronoi mesh across MPI ranks by graph
partition using METIS {cite:p}`karypis1998metis`, with each rank
owning a contiguous sub-graph of cells selected to balance load
and minimize edge cut. SHDOM, whose internal algorithm is
naturally Cartesian, decomposes its domain differently: it
establishes a 2D $(X, Y)$ pencil decomposition via
`MPI_DIMS_CREATE` and `MPI_CART_CREATE`, with the vertical
dimension $Z$ kept local on each rank. The two partitions have no
geometric relation to each other, so optical-property data must be
redistributed from MPAS-A's graph partition onto SHDOM's pencil
partition before each SHDOM call, and the returned fluxes must be
redistributed back.

The framework executes this redistribution through a precomputed
`MPI_Alltoallv` plan. At framework initialization each rank
mirrors the `MPI_DIMS_CREATE` logic of SHDOM on the MPAS-A side,
arriving at the same pencil tiling; the mapper's sparse weight
structure is then augmented with, for each Cartesian cell, the
rank index on which that cell lives in the SHDOM decomposition.
From these per-cell rank assignments the collective `MPI_Alltoallv`
counts and displacements are computed once and retained. At each
radiation step in which SHDOM runs, one pair of `MPI_Alltoallv`
calls implements the forward and reverse redistributions. For the
Phase 1 supercell configuration the payload is a few megabytes per
radiation step (tens of thousands of cells, tens of levels, three
optical-property components, a single band), and the
redistribution is not a performance bottleneck.

A small compatibility patch is required on the SHDOM side.
Upstream SHDOM was designed as a standalone executable and
hard-codes its communicator as `MPI_COMM_WORLD` at fourteen call
sites, and calls `MPI_Init` unconditionally during its setup
phase. Running inside the MPAS-A process, which has already
initialized MPI and which may wish to use a communicator other
than the world, those assumptions must be relaxed. The patch —
some fifteen lines, documented in
[](../appendices/shdom-mpi-patch.md) — adds an optional
communicator argument to `START_MPI` and `MAP_SHDOM_MPI`, skips
`MPI_Init` when a communicator is supplied, and replaces the
fourteen hard-coded `MPI_COMM_WORLD` references with a
module-scope variable populated from the supplied communicator.
Standalone SHDOM continues to function unchanged: if the patched
routines are called with no communicator argument they recover the
original behaviour. The patch is maintained as a tracked commit on
the vendored SHDOM subtree.
