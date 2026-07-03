program test_elliptic_operator
   use mpas_kind_types, only : RKIND
   use mpas_elliptic_operator, only : mpas_elliptic_operator_type, elliptic_operator_init, &
                                      elliptic_operator_init_from_weights
   implicit none

   call test_constant_has_ground_dirichlet_only()
   call test_diagonal_is_positive()
   call test_horizontal_loop_uses_nEdgesOnCell()
   call test_init_from_weights_identity_walls_symmetry()
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

   subroutine test_init_from_weights_identity_walls_symmetry()
      type(mpas_elliptic_operator_type) :: op
      integer, parameter :: nCells = 2, nEdges = 1, nVertLevels = 3, maxEdges = 1
      integer :: nEdgesOnCell(nCells), cellsOnCell(maxEdges,nCells), edgesOnCell(maxEdges,nCells)
      real(kind=RKIND) :: h_weight(nVertLevels,nEdges)
      real(kind=RKIND) :: v_up(nVertLevels,nCells), v_lo(nVertLevels,nCells)
      real(kind=RKIND) :: diag_extra(nVertLevels,nCells), volume(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)
      real(kind=RKIND) :: x(nVertLevels,nCells), y(nVertLevels,nCells), y12, y21

      nEdgesOnCell = [1, 1]
      cellsOnCell(:,1) = [2]
      cellsOnCell(:,2) = [1]
      edgesOnCell(:,1) = [1]
      edgesOnCell(:,2) = [1]

      h_weight(:,1) = [0.0_RKIND, 0.8_RKIND, 1.0_RKIND]
      v_up(:,1) = [1.0_RKIND, 1.0_RKIND, 0.0_RKIND]
      v_lo(:,1) = [2.0_RKIND, 1.0_RKIND, 1.0_RKIND]
      diag_extra(:,1) = [2.0_RKIND, 0.4_RKIND, 0.0_RKIND]
      volume(:,1) = [1.0_RKIND, 1.0_RKIND, 1.0_RKIND]
      active(:,1) = .true.
      v_up(:,2) = [0.0_RKIND, 1.0_RKIND/0.9_RKIND, 0.0_RKIND]
      v_lo(:,2) = [0.0_RKIND, 2.5_RKIND, 1.0_RKIND/0.9_RKIND]
      diag_extra(:,2) = 0.0_RKIND
      volume(:,2) = [0.0_RKIND, 0.8_RKIND, 1.0_RKIND]
      active(:,2) = [.false., .true., .true.]

      call elliptic_operator_init_from_weights(op, nCells, nEdges, nVertLevels, maxEdges, &
                                               nEdgesOnCell, cellsOnCell, edgesOnCell, &
                                               h_weight, v_up, v_lo, diag_extra, volume, active)

      x = 0.0_RKIND
      x(1,2) = 1.0_RKIND
      call op%matvec(x, y)
      if (abs(y(1,2) - 1.0_RKIND) > 1.0e-14_RKIND) stop "FAIL: inactive row not identity"
      if (any(abs(y(:,1)) > 1.0e-14_RKIND) .or. any(abs(y(2:3,2)) > 1.0e-14_RKIND)) &
         stop "FAIL: inactive cell leaks coupling"

      x = 0.0_RKIND
      x(1,1) = 1.0_RKIND
      call op%matvec(x, y)
      if (abs(y(1,1) - 5.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: diag with wall sink"
      if (abs(y(2,1) + 1.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: vertical coupling"

      x = 0.0_RKIND
      x(2,2) = 1.0_RKIND
      call op%matvec(x, y)
      y21 = y(2,1)
      x = 0.0_RKIND
      x(2,1) = 1.0_RKIND
      call op%matvec(x, y)
      y12 = y(2,2)
      if (abs(y21 + 0.8_RKIND) > 1.0e-12_RKIND .or. abs(y21 - y12) > 1.0e-15_RKIND) &
         stop "FAIL: open-face coupling not symmetric"

      x = 0.0_RKIND
      where (active) x = 1.0_RKIND
      call op%matvec(x, y)
      if (abs(y(2,2) - 2.5_RKIND) > 1.0e-12_RKIND) stop "FAIL: cut-cell ground sink"
      if (any(abs(y(3,:)) > 1.0e-12_RKIND)) stop "FAIL: interior rows should annihilate constants"

      if (abs(op%volume(2,2) - 0.8_RKIND) > 1.0e-14_RKIND) stop "FAIL: volume not copied"
      call op%destroy()
   end subroutine test_init_from_weights_identity_walls_symmetry

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
