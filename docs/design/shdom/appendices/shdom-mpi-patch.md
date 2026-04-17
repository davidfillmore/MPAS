(app-shdom-mpi-patch)=
# SHDOM MPI patch

SHDOM was designed as a standalone executable and calls
`MPI_Init` unconditionally from its setup routine `START_MPI` in
`shdom_mpi.f90`. A further fourteen call sites in the same file
reference `MPI_COMM_WORLD` directly. Embedding SHDOM inside the
MPAS-A process, in which MPI has already been initialized and in
which the working communicator may be a subset rather than the
world communicator, requires that both assumptions be relaxed
without disturbing SHDOM's standalone behaviour.

The patch is approximately fifteen lines, applied as a tracked
commit on the vendored SHDOM subtree. Three changes suffice:

1. An optional `comm_in` argument is added to `START_MPI` and to
   `MAP_SHDOM_MPI`.
2. When `comm_in` is present and MPI has already been initialized,
   `MPI_Init` is skipped and the supplied communicator is stored
   in a module-scope variable `comm2d`.
3. The fourteen direct references to `MPI_COMM_WORLD` are replaced
   by `comm2d`, which defaults to `MPI_COMM_WORLD` when no
   external communicator is supplied. Standalone SHDOM invocations,
   which do not pass `comm_in`, recover the original behaviour
   exactly.

The fourteen call sites are audited by searching for
`MPI_COMM_WORLD` in `shdom_mpi.f90` and confirming that each
occurrence has been replaced by `comm2d`; this audit is the
correctness proof for the backward-compatibility claim. The
patch is reapplied as a rebase on any future SHDOM upstream
update.
