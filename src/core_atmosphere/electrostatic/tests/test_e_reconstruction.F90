program test_e_reconstruction

   use mpas_kind_types, only : RKIND
   use mpas_electrostatic_bcs, only : electrostatic_apply_ground_rhs, &
                                      electrostatic_compute_vertical_E, &
                                      electrostatic_horizontal_metric_scale

   implicit none

   call test_vertical_e_ground_second_order_and_top_neumann()
   call test_ground_rhs_only_modifies_bottom_level()
   call test_horizontal_metric_scale_respects_mesh_geometry()
   print *, "PASS: E reconstruction tests"

contains

   subroutine test_vertical_e_ground_second_order_and_top_neumann()
      integer, parameter :: nCells = 1, nVertLevels = 3
      real(kind=RKIND) :: phi(nVertLevels,nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: E_vertical(nVertLevels,nCells)
      real(kind=RKIND) :: expected(nVertLevels,nCells), phi_ground, z_top
      real(kind=RKIND) :: zGround(nCells)
      integer :: kFirst(nCells)

      ! phi = z^2 sampled at the layer midpoints, consistent with phi=0 at the
      ! ground face (z=0). The ground-adjacent estimate is a second-order
      ! one-sided derivative, exact for a quadratic: d/dz(z^2) = 2z.
      zMid(:,1) = [0.5_RKIND, 1.5_RKIND, 2.5_RKIND]
      phi(:,1) = zMid(:,1)**2
      phi_ground = 0.0_RKIND
      z_top = 3.0_RKIND
      zGround = 0.0_RKIND
      kFirst = 1

      call electrostatic_compute_vertical_E(nCells, nVertLevels, phi, zMid, phi_ground, z_top, &
                                            zGround, kFirst, E_vertical)

      expected(1,1) = -2.0_RKIND * zMid(1,1)
      expected(2,1) = -(phi(3,1) - phi(1,1)) / (zMid(3,1) - zMid(1,1))
      expected(3,1) = -2.0_RKIND * (phi(2,1) - phi(3,1)) * (zMid(3,1) - z_top) / &
                      ((zMid(2,1) - z_top)**2 - (zMid(3,1) - z_top)**2)

      if (maxval(abs(E_vertical - expected)) > 1.0e-12_RKIND) stop "FAIL: vertical E reconstruction mismatch"

      zGround = 1.0_RKIND
      kFirst = 2
      phi(:,1) = 0.0_RKIND
      phi(2,1) = (zMid(2,1) - 1.0_RKIND)**2
      phi(3,1) = (zMid(3,1) - 1.0_RKIND)**2
      call electrostatic_compute_vertical_E(nCells, nVertLevels, phi, zMid, phi_ground, z_top, &
                                            zGround, kFirst, E_vertical)
      if (abs(E_vertical(2,1) + 2.0_RKIND*(zMid(2,1)-1.0_RKIND)) > 1.0e-12_RKIND) &
         stop "FAIL: lifted-ground quadratic"
      if (abs(E_vertical(1,1)) > 0.0_RKIND) stop "FAIL: buried level E /= 0"
   end subroutine test_vertical_e_ground_second_order_and_top_neumann

   subroutine test_ground_rhs_only_modifies_bottom_level()
      integer, parameter :: nCells = 2, nVertLevels = 3
      real(kind=RKIND) :: rhs(nVertLevels,nCells), v_weight_lower(nVertLevels,nCells)

      rhs = 0.0_RKIND
      v_weight_lower = 0.0_RKIND
      v_weight_lower(1,:) = [2.0_RKIND, 3.0_RKIND]

      call electrostatic_apply_ground_rhs(nCells, nVertLevels, 5.0_RKIND, v_weight_lower, rhs)

      if (abs(rhs(1,1) - 10.0_RKIND) > 1.0e-14_RKIND) stop "FAIL: ground RHS first cell mismatch"
      if (abs(rhs(1,2) - 15.0_RKIND) > 1.0e-14_RKIND) stop "FAIL: ground RHS second cell mismatch"
      if (any(abs(rhs(2:nVertLevels,:)) > 1.0e-14_RKIND)) stop "FAIL: ground RHS changed non-ground levels"
   end subroutine test_ground_rhs_only_modifies_bottom_level

   subroutine test_horizontal_metric_scale_respects_mesh_geometry()
      real(kind=RKIND), parameter :: radius = 6371229.0_RKIND
      real(kind=RKIND) :: planar_dc(2), unit_sphere_dc(2), meter_sphere_dc(2)

      planar_dc = [500.0_RKIND, 750.0_RKIND]
      unit_sphere_dc = [0.01_RKIND, 0.02_RKIND]
      meter_sphere_dc = [50000.0_RKIND, 75000.0_RKIND]

      if (electrostatic_horizontal_metric_scale(2, planar_dc, radius, .false.) /= 1.0_RKIND) &
         stop "FAIL: planar mesh metrics should not be radius-scaled"

      if (electrostatic_horizontal_metric_scale(2, unit_sphere_dc, radius, .true.) /= radius) &
         stop "FAIL: unit-sphere mesh metrics should be radius-scaled"

      if (electrostatic_horizontal_metric_scale(2, meter_sphere_dc, radius, .true.) /= 1.0_RKIND) &
         stop "FAIL: meter-scale spherical mesh metrics should not be radius-scaled"
   end subroutine test_horizontal_metric_scale_respects_mesh_geometry

end program test_e_reconstruction
