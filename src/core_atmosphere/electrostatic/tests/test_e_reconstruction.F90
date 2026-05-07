program test_e_reconstruction

   use mpas_kind_types, only : RKIND
   use mpas_electrostatic_bcs, only : electrostatic_apply_ground_rhs, &
                                      electrostatic_compute_vertical_E

   implicit none

   call test_vertical_e_uses_ground_centered_and_top_neumann()
   call test_ground_rhs_only_modifies_bottom_level()
   print *, "PASS: E reconstruction tests"

contains

   subroutine test_vertical_e_uses_ground_centered_and_top_neumann()
      integer, parameter :: nCells = 1, nVertLevels = 3
      real(kind=RKIND) :: phi(nVertLevels,nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: E_vertical(nVertLevels,nCells)
      real(kind=RKIND) :: expected(nVertLevels,nCells), phi_ground

      phi(:,1) = [0.0_RKIND, 2.0_RKIND, 4.0_RKIND]
      zMid(:,1) = [0.5_RKIND, 1.5_RKIND, 2.5_RKIND]
      phi_ground = 0.0_RKIND

      call electrostatic_compute_vertical_E(nCells, nVertLevels, phi, zMid, phi_ground, E_vertical)

      expected(1,1) = -(phi(2,1) - phi_ground) / (zMid(2,1) - 0.0_RKIND)
      expected(2,1) = -(phi(3,1) - phi(1,1)) / (zMid(3,1) - zMid(1,1))
      expected(3,1) = 0.0_RKIND

      if (maxval(abs(E_vertical - expected)) > 1.0e-14_RKIND) stop "FAIL: vertical E reconstruction mismatch"
   end subroutine test_vertical_e_uses_ground_centered_and_top_neumann

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

end program test_e_reconstruction
