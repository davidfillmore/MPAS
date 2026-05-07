program test_pcg_solver

   use mpas_kind_types, only : RKIND
   use mpas_elliptic_operator, only : mpas_elliptic_operator_type, elliptic_operator_init
   use mpas_pcg_solver, only : mpas_pcg_solver_type, pcg_init, PC_JACOBI, PC_NONE

   implicit none

   call test_zero_rhs_converges_without_iterations()
   call test_recovers_manufactured_solution()
   print *, "PASS: PCG solver tests"

contains

   subroutine test_zero_rhs_converges_without_iterations()
      type(mpas_elliptic_operator_type) :: op
      type(mpas_pcg_solver_type) :: solver
      real(kind=RKIND) :: rhs(3,2), phi(3,2), r0, rfinal
      integer :: iter_count, status

      call make_two_cell_operator(op)
      call pcg_init(solver, 2, 3, PC_NONE, 5, 1.0e-12_RKIND)

      rhs = 0.0_RKIND
      phi = 0.0_RKIND
      call solver%solve(op, rhs, phi, iter_count, r0, rfinal, status)

      if (status /= 0) stop "FAIL: zero RHS did not converge"
      if (iter_count /= 0) stop "FAIL: zero RHS should require no PCG iterations"
      if (r0 > 1.0e-24_RKIND) stop "FAIL: zero RHS initial residual is nonzero"
      if (rfinal > 1.0e-24_RKIND) stop "FAIL: zero RHS final residual is nonzero"
      if (any(abs(phi) > 1.0e-24_RKIND)) stop "FAIL: zero RHS changed the solution"

      call solver%destroy()
      call op%destroy()
   end subroutine test_zero_rhs_converges_without_iterations

   subroutine test_recovers_manufactured_solution()
      type(mpas_elliptic_operator_type) :: op
      type(mpas_pcg_solver_type) :: solver
      real(kind=RKIND) :: rhs(3,2), phi(3,2), phi_exact(3,2), r0, rfinal
      integer :: iter_count, status

      call make_two_cell_operator(op)
      call pcg_init(solver, 2, 3, PC_JACOBI, 20, 1.0e-12_RKIND)

      phi_exact(:,1) = [0.125_RKIND, 0.25_RKIND, 0.5_RKIND]
      phi_exact(:,2) = [-0.25_RKIND, 0.375_RKIND, 0.625_RKIND]
      call op%matvec(phi_exact, rhs)

      phi = 0.0_RKIND
      call solver%solve(op, rhs, phi, iter_count, r0, rfinal, status)

      if (status /= 0) stop "FAIL: manufactured solve did not converge"
      if (iter_count <= 0 .or. iter_count > 20) stop "FAIL: unexpected iteration count"
      if (r0 <= 0.0_RKIND) stop "FAIL: manufactured solve initial residual not positive"
      if (rfinal > 1.0e-10_RKIND) stop "FAIL: manufactured solve residual too large"
      if (maxval(abs(phi - phi_exact)) > 1.0e-8_RKIND) stop "FAIL: manufactured solution mismatch"

      call solver%destroy()
      call op%destroy()
   end subroutine test_recovers_manufactured_solution

   subroutine make_two_cell_operator(op)
      type(mpas_elliptic_operator_type), intent(out) :: op
      integer, parameter :: nCells = 2, nEdges = 1, nVertLevels = 3, maxEdges = 1
      integer :: nEdgesOnCell(nCells)
      integer :: cellsOnCell(maxEdges,nCells), edgesOnCell(maxEdges,nCells)
      real(kind=RKIND) :: areaCell(nCells), dcEdge(nEdges), dvEdge(nEdges), dz(nVertLevels)

      nEdgesOnCell = [1, 1]
      cellsOnCell(:,1) = [2]
      cellsOnCell(:,2) = [1]
      edgesOnCell(:,1) = [1]
      edgesOnCell(:,2) = [1]
      areaCell = 1.0_RKIND
      dcEdge = 1.0_RKIND
      dvEdge = 1.0_RKIND
      dz = 1.0_RKIND

      call elliptic_operator_init(op, nCells, nEdges, nVertLevels, maxEdges, &
                                  areaCell, dcEdge, dvEdge, nEdgesOnCell, &
                                  cellsOnCell, edgesOnCell, dz, 8.854e-12_RKIND)
   end subroutine make_two_cell_operator

end program test_pcg_solver
