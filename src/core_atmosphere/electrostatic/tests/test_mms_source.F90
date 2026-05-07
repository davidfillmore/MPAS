program test_mms_source

   use mpas_kind_types, only : RKIND
   use mpas_electrostatic_benchmarks, only : electrostatic_mms_cart_phi_exact, &
                                             electrostatic_compute_error_norms
   use mpas_electrostatic_source, only : electrostatic_fill_mms_cart_source

   implicit none

   call test_mms_source_uses_sine_vertical_basis()
   call test_error_norms_are_relative_l2_and_absolute_linf()
   print *, "PASS: MMS source utility tests"

contains

   subroutine test_mms_source_uses_sine_vertical_basis()
      integer, parameter :: nCells = 2, nVertLevels = 2
      real(kind=RKIND), parameter :: epsilon0 = 8.8541878128e-12_RKIND
      real(kind=RKIND), parameter :: x_period = 2.0_RKIND
      real(kind=RKIND), parameter :: y_period = 4.0_RKIND
      real(kind=RKIND), parameter :: z_top = 1.0_RKIND
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND) :: phi_exact, coeff, expected, top_deriv
      real(kind=RKIND) :: kx, ky, kz, h

      xCell = [0.5_RKIND, 1.0_RKIND]
      yCell = [1.0_RKIND, 2.0_RKIND]
      zMid(:,1) = [0.25_RKIND, 0.75_RKIND]
      zMid(:,2) = [0.25_RKIND, 0.75_RKIND]

      call electrostatic_fill_mms_cart_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                              epsilon0, x_period, y_period, z_top, rho_charge)

      kx = 2.0_RKIND * acos(-1.0_RKIND) / x_period
      ky = 2.0_RKIND * acos(-1.0_RKIND) / y_period
      kz = 0.5_RKIND * acos(-1.0_RKIND) / z_top
      coeff = epsilon0 * (kx*kx + ky*ky + kz*kz)
      phi_exact = electrostatic_mms_cart_phi_exact(xCell(1), yCell(1), zMid(1,1), &
                                                   x_period, y_period, z_top)
      expected = coeff * phi_exact

      if (abs(rho_charge(1,1) - expected) > 1.0e-20_RKIND) stop "FAIL: MMS charge source mismatch"
      if (abs(electrostatic_mms_cart_phi_exact(xCell(1), yCell(1), 0.0_RKIND, &
                                               x_period, y_period, z_top)) > 1.0e-14_RKIND) &
         stop "FAIL: MMS exact potential is not grounded"

      h = 1.0e-6_RKIND
      top_deriv = (electrostatic_mms_cart_phi_exact(xCell(1), yCell(1), z_top + h, &
                                                    x_period, y_period, z_top) - &
                   electrostatic_mms_cart_phi_exact(xCell(1), yCell(1), z_top - h, &
                                                    x_period, y_period, z_top)) / (2.0_RKIND * h)
      if (abs(top_deriv) > 1.0e-9_RKIND) stop "FAIL: MMS exact potential violates top Neumann condition"
   end subroutine test_mms_source_uses_sine_vertical_basis

   subroutine test_error_norms_are_relative_l2_and_absolute_linf()
      real(kind=RKIND) :: phi_exact(2,2), phi(2,2), volume(2,2)
      real(kind=RKIND) :: l2, linf

      phi_exact(:,1) = [1.0_RKIND, 2.0_RKIND]
      phi_exact(:,2) = [3.0_RKIND, 4.0_RKIND]
      phi = 2.0_RKIND * phi_exact
      volume = 1.0_RKIND

      call electrostatic_compute_error_norms(phi, phi_exact, volume, l2, linf)

      if (abs(l2 - 1.0_RKIND) > 1.0e-14_RKIND) stop "FAIL: MMS relative L2 norm mismatch"
      if (abs(linf - 4.0_RKIND) > 1.0e-14_RKIND) stop "FAIL: MMS Linf norm mismatch"
   end subroutine test_error_norms_are_relative_l2_and_absolute_linf

end program test_mms_source
