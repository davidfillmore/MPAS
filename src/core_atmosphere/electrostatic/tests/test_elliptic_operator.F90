program test_elliptic_operator
   use mpas_kind_types, only : RKIND
   use mpas_elliptic_operator, only : mpas_elliptic_operator_type, elliptic_operator_init
   implicit none

   call test_constant_has_ground_dirichlet_only()
   call test_diagonal_is_positive()
   call test_horizontal_loop_uses_nEdgesOnCell()
   print *, "PASS: elliptic operator tests"

contains

   subroutine test_constant_has_ground_dirichlet_only()
      type(mpas_elliptic_operator_type) :: op
      real(kind=RKIND) :: x(3,2), y(3,2)

      call make_two_cell_operator(op)
      x = 1.0_RKIND
      call op%matvec(x, y)

      if (any(abs(y(2:3,:)) > 1.0e-24_RKIND)) stop "FAIL: constant changed away from ground row"
      if (.not. all(y(1,:) > 0.0_RKIND)) stop "FAIL: ground Dirichlet not in A diagonal"

      call op%destroy()
   end subroutine test_constant_has_ground_dirichlet_only

   subroutine test_diagonal_is_positive()
      type(mpas_elliptic_operator_type) :: op

      call make_two_cell_operator(op)
      if (.not. all(op%diag > 0.0_RKIND)) stop "FAIL: A diagonal must be positive"

      call op%destroy()
   end subroutine test_diagonal_is_positive

   subroutine test_horizontal_loop_uses_nEdgesOnCell()
      type(mpas_elliptic_operator_type) :: op
      real(kind=RKIND) :: x(3,2), y(3,2)

      call make_two_cell_operator(op)
      op%nEdgesOnCell(1) = 1
      op%nEdgesOnCell(2) = 0
      x = 0.0_RKIND
      x(:,1) = 1.0_RKIND

      call op%matvec(x, y)
      if (any(abs(y(:,2)) > 1.0e-24_RKIND)) stop "FAIL: matvec ignored nEdgesOnCell"

      call op%destroy()
   end subroutine test_horizontal_loop_uses_nEdgesOnCell

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

end program test_elliptic_operator
