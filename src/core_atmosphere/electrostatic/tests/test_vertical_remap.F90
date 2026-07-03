program test_vertical_remap

   use mpas_kind_types, only : RKIND
   use mpas_electrostatic_vertical_remap, only : electrostatic_interp_mpas_to_poisson, &
                                                electrostatic_interp_w_to_poisson_mid, &
                                                electrostatic_build_poisson_grid, &
                                                electrostatic_build_poisson_fv_weights

   implicit none

   call test_linear_profile_exact()
   call test_inactive_poisson_levels_zero()
   call test_w_interfaces_interpolate_to_poisson_midpoints()
   call test_build_poisson_grid_flat_matches_current_operator()
   call test_build_poisson_grid_terrain_clipping_and_sliver_merge()
   call test_build_poisson_fv_weights_shaved_and_grounded()
   print *, "PASS: vertical remap tests"

contains

   subroutine test_linear_profile_exact()
      integer, parameter :: nCells = 1, nVertLevels = 4
      real(kind=RKIND) :: zMpas(nVertLevels,nCells), zPoisson(nVertLevels,nCells)
      real(kind=RKIND) :: source(nVertLevels,nCells), target(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      zMpas(:,1) = [100.0_RKIND, 300.0_RKIND, 700.0_RKIND, 900.0_RKIND]
      source(:,1) = 2.0_RKIND * zMpas(:,1) + 5.0_RKIND
      zPoisson(:,1) = [100.0_RKIND, 200.0_RKIND, 500.0_RKIND, 900.0_RKIND]
      active(:,1) = .true.

      call electrostatic_interp_mpas_to_poisson(nCells, nVertLevels, zMpas, zPoisson, active, source, target)

      if (maxval(abs(target(:,1) - (2.0_RKIND*zPoisson(:,1) + 5.0_RKIND))) > 1.0e-12_RKIND) &
         stop "FAIL: linear interpolation not exact"
   end subroutine test_linear_profile_exact

   subroutine test_inactive_poisson_levels_zero()
      integer, parameter :: nCells = 1, nVertLevels = 3
      real(kind=RKIND) :: zMpas(nVertLevels,nCells), zPoisson(nVertLevels,nCells)
      real(kind=RKIND) :: source(nVertLevels,nCells), target(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      zMpas(:,1) = [100.0_RKIND, 300.0_RKIND, 500.0_RKIND]
      source(:,1) = [10.0_RKIND, 30.0_RKIND, 50.0_RKIND]
      zPoisson(:,1) = [50.0_RKIND, 100.0_RKIND, 200.0_RKIND]
      active(:,1) = [.false., .true., .true.]

      call electrostatic_interp_mpas_to_poisson(nCells, nVertLevels, zMpas, zPoisson, active, source, target)

      if (target(1,1) /= 0.0_RKIND) stop "FAIL: inactive level not zeroed"
      if (abs(target(2,1) - 10.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: lower clamp mismatch"
      if (abs(target(3,1) - 20.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: interior interpolation mismatch"
   end subroutine test_inactive_poisson_levels_zero

   subroutine test_w_interfaces_interpolate_to_poisson_midpoints()
      integer, parameter :: nCells = 1, nVertLevels = 2
      real(kind=RKIND) :: zgrid(nVertLevels+1,nCells), zPoisson(nVertLevels,nCells)
      real(kind=RKIND) :: w(nVertLevels+1,nCells), wPoisson(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      zgrid(:,1) = [0.0_RKIND, 100.0_RKIND, 300.0_RKIND]
      w(:,1) = [0.0_RKIND, 10.0_RKIND, 30.0_RKIND]
      zPoisson(:,1) = [50.0_RKIND, 200.0_RKIND]
      active(:,1) = .true.

      call electrostatic_interp_w_to_poisson_mid(nCells, nVertLevels, zgrid, zPoisson, active, w, wPoisson)

      if (abs(wPoisson(1,1) - 5.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: w interpolation level 1"
      if (abs(wPoisson(2,1) - 20.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: w interpolation level 2"
   end subroutine test_w_interfaces_interpolate_to_poisson_midpoints

   subroutine test_build_poisson_grid_flat_matches_current_operator()
      integer, parameter :: nCells = 1, nVertLevels = 4
      real(kind=RKIND) :: zGround(nCells), zFace(nVertLevels+1)
      real(kind=RKIND) :: poissonZmid(nVertLevels,nCells), airThickness(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)
      integer :: kFirst(nCells)

      zGround = 0.0_RKIND
      call electrostatic_build_poisson_grid(nCells, nVertLevels, zGround, 0.0_RKIND, 8.0_RKIND, &
                                            zFace, poissonZmid, active, kFirst, airThickness)
      if (abs(zFace(1)) > 0.0_RKIND .or. abs(zFace(nVertLevels+1) - 8.0_RKIND) > 1.0e-12_RKIND) &
         stop "FAIL: face range"
      if (kFirst(1) /= 1) stop "FAIL: flat column should start at level 1"
      if (.not. all(active(:,1))) stop "FAIL: flat column fully active"
      if (maxval(abs(airThickness(:,1) - 2.0_RKIND)) > 1.0e-12_RKIND) stop "FAIL: flat air thickness"
      if (maxval(abs(poissonZmid(:,1) - [1.0_RKIND, 3.0_RKIND, 5.0_RKIND, 7.0_RKIND])) > &
          1.0e-12_RKIND) stop "FAIL: flat midpoints must be uniform band midpoints"
   end subroutine test_build_poisson_grid_flat_matches_current_operator

   subroutine test_build_poisson_grid_terrain_clipping_and_sliver_merge()
      integer, parameter :: nCells = 2, nVertLevels = 4
      real(kind=RKIND) :: zGround(nCells), zFace(nVertLevels+1)
      real(kind=RKIND) :: poissonZmid(nVertLevels,nCells), airThickness(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)
      integer :: kFirst(nCells)

      ! dzP = 2: faces 0,2,4,6,8
      ! cell 1: zg=2.5 -> band 2 air = 1.5, kFirst=2, centroid 3.25.
      ! cell 2: zg=3.95 -> sliver band 2 merges up, kFirst=3, centroid 4.975.
      zGround = [2.5_RKIND, 3.95_RKIND]
      call electrostatic_build_poisson_grid(nCells, nVertLevels, zGround, 0.0_RKIND, 8.0_RKIND, &
                                            zFace, poissonZmid, active, kFirst, airThickness)
      if (kFirst(1) /= 2 .or. kFirst(2) /= 3) stop "FAIL: kFirstActive"
      if (active(1,1) .or. active(2,2)) stop "FAIL: below-terrain levels must be inactive"
      if (abs(airThickness(2,1) - 1.5_RKIND) > 1.0e-12_RKIND) stop "FAIL: kFirst band cell 1"
      if (abs(poissonZmid(2,1) - 3.25_RKIND) > 1.0e-12_RKIND) stop "FAIL: air-centroid midpoint cell 1"
      if (abs(airThickness(3,2) - 2.05_RKIND) > 1.0e-12_RKIND) stop "FAIL: merged band cell 2"
      if (abs(poissonZmid(3,2) - 4.975_RKIND) > 1.0e-12_RKIND) stop "FAIL: merged centroid cell 2"
      if (abs(airThickness(1,1)) > 0.0_RKIND) stop "FAIL: submerged level thickness"
      if (abs(sum(airThickness(:,1)) - (8.0_RKIND - 2.5_RKIND)) > 1.0e-12_RKIND) &
         stop "FAIL: air-column sum cell 1"
      if (abs(sum(airThickness(:,2)) - (8.0_RKIND - 3.95_RKIND)) > 1.0e-12_RKIND) &
         stop "FAIL: air-column sum cell 2"
   end subroutine test_build_poisson_grid_terrain_clipping_and_sliver_merge

   subroutine test_build_poisson_fv_weights_shaved_and_grounded()
      integer, parameter :: nCells = 2, nEdges = 1, nVertLevels = 3
      integer :: cellsOnEdge(2,nEdges), kFirst(nCells)
      real(kind=RKIND) :: areaCell(nCells), dcEdge(nEdges), dvEdge(nEdges)
      real(kind=RKIND) :: zGround(nCells), zFace(nVertLevels+1)
      real(kind=RKIND) :: poissonZmid(nVertLevels,nCells), airThickness(nVertLevels,nCells)
      real(kind=RKIND) :: h_weight(nVertLevels,nEdges)
      real(kind=RKIND) :: v_up(nVertLevels,nCells), v_lo(nVertLevels,nCells)
      real(kind=RKIND) :: diag_extra(nVertLevels,nCells), volume(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      cellsOnEdge(:,1) = [1, 2]
      areaCell = 1.0_RKIND
      dcEdge = 1.0_RKIND
      dvEdge = 1.0_RKIND
      zGround = [0.0_RKIND, 1.2_RKIND]
      call electrostatic_build_poisson_grid(nCells, nVertLevels, zGround, 0.0_RKIND, 3.0_RKIND, &
                                            zFace, poissonZmid, active, kFirst, airThickness)
      call electrostatic_build_poisson_fv_weights(nCells, nEdges, nVertLevels, cellsOnEdge, &
                                                  areaCell, dcEdge, dvEdge, 1.0_RKIND, &
                                                  zFace, zGround, 3.0_RKIND, poissonZmid, &
                                                  active, kFirst, airThickness, &
                                                  h_weight, v_up, v_lo, diag_extra, volume)

      if (abs(h_weight(1,1)) > 0.0_RKIND) stop "FAIL: band 1 must be closed"
      if (abs(h_weight(2,1) - 0.8_RKIND) > 1.0e-12_RKIND) stop "FAIL: band 2 open weight"
      if (abs(h_weight(3,1) - 1.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: band 3 open weight"
      if (abs(diag_extra(1,1) - 2.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: band 1 wall sink"
      if (abs(diag_extra(2,1) - 0.4_RKIND) > 1.0e-12_RKIND) stop "FAIL: band 2 wall sink"
      if (any(abs(diag_extra(:,2)) > 0.0_RKIND)) stop "FAIL: higher-ground cell has no wall"
      if (abs(v_lo(1,1) - 2.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: flat ground weight"
      if (abs(v_lo(2,2) - 2.5_RKIND) > 1.0e-12_RKIND) stop "FAIL: cut-cell ground weight"
      if (abs(v_up(2,2) - 1.0_RKIND/0.9_RKIND) > 1.0e-12_RKIND) stop "FAIL: centroid vertical weight"
      if (abs(v_lo(3,2) - v_up(2,2)) > 1.0e-15_RKIND) stop "FAIL: vertical weight pair symmetry"
      if (abs(volume(2,2) - 0.8_RKIND) > 1.0e-12_RKIND) stop "FAIL: volume /= area*airThickness"
      if (abs(sum(h_weight(:,1)) + 0.5_RKIND*sum(diag_extra(:,1)) - 3.0_RKIND) > 1.0e-12_RKIND) &
         stop "FAIL: face partition does not tile the air column"
   end subroutine test_build_poisson_fv_weights_shaved_and_grounded

end program test_vertical_remap
