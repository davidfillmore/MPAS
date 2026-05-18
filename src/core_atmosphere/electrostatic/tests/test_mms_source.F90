program test_mms_source

   use mpas_kind_types, only : RKIND
   use mpas_electrostatic_benchmarks, only : electrostatic_mms_cart_phi_exact, &
                                             electrostatic_compute_error_norms
   use mpas_electrostatic_source, only : electrostatic_fill_mms_cart_source, &
                                         electrostatic_fill_mms_cart_horizontal_source, &
                                         electrostatic_fill_mms_sphere_source, &
                                         electrostatic_fill_mms_sphere_horizontal_source, &
                                         electrostatic_fill_dipole_source, &
                                         electrostatic_fill_tripole_source, &
                                         electrostatic_fill_gaussian_monopole_source, &
                                         electrostatic_fill_gaussian_dipole_y_source, &
                                         electrostatic_fill_gaussian_dipole_z_source

   implicit none

   call test_mms_source_uses_sine_vertical_basis()
   call test_horizontal_mms_source_uses_discrete_vertical_operator()
   call test_sphere_mms_source_uses_boundary_compatible_vertical_basis()
   call test_sphere_horizontal_mms_source_uses_discrete_vertical_operator()
   call test_dipole_source_uses_explicit_center_with_expected_lobe_signs()
   call test_tripole_source_uses_explicit_center_with_expected_lobe_signs()
   call test_gaussian_monopole_source_is_centered_and_positive()
   call test_gaussian_dipole_sources_have_expected_signs()
   call test_gaussian_charges_have_expected_integrals()
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

   subroutine test_horizontal_mms_source_uses_discrete_vertical_operator()
      integer, parameter :: nCells = 1, nVertLevels = 2
      real(kind=RKIND), parameter :: epsilon0 = 8.8541878128e-12_RKIND
      real(kind=RKIND), parameter :: x_period = 2.0_RKIND
      real(kind=RKIND), parameter :: y_period = 4.0_RKIND
      real(kind=RKIND), parameter :: z_top = 1.0_RKIND
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND) :: v_weight_upper(nVertLevels,nCells)
      real(kind=RKIND) :: v_weight_lower(nVertLevels,nCells)
      real(kind=RKIND) :: volume(nVertLevels,nCells)
      real(kind=RKIND) :: phi1, phi2, horizontal_coeff, expected1, expected2
      real(kind=RKIND) :: kx, ky

      xCell = [0.5_RKIND]
      yCell = [1.0_RKIND]
      zMid(:,1) = [0.25_RKIND, 0.75_RKIND]
      volume(:,1) = [2.0_RKIND, 3.0_RKIND]
      v_weight_upper(:,1) = [5.0_RKIND, 0.0_RKIND]
      v_weight_lower(:,1) = [7.0_RKIND, 11.0_RKIND]

      call electrostatic_fill_mms_cart_horizontal_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                         epsilon0, x_period, y_period, z_top, &
                                                         v_weight_upper, v_weight_lower, volume, &
                                                         rho_charge)

      kx = 2.0_RKIND * acos(-1.0_RKIND) / x_period
      ky = 2.0_RKIND * acos(-1.0_RKIND) / y_period
      horizontal_coeff = epsilon0 * (kx*kx + ky*ky)
      phi1 = electrostatic_mms_cart_phi_exact(xCell(1), yCell(1), zMid(1,1), &
                                              x_period, y_period, z_top)
      phi2 = electrostatic_mms_cart_phi_exact(xCell(1), yCell(1), zMid(2,1), &
                                              x_period, y_period, z_top)

      expected1 = horizontal_coeff * phi1 + &
                  (v_weight_upper(1,1) * (phi1 - phi2) + &
                   v_weight_lower(1,1) * phi1) / volume(1,1)
      expected2 = horizontal_coeff * phi2 + &
                  (v_weight_lower(2,1) * (phi2 - phi1)) / volume(2,1)

      if (abs(rho_charge(1,1) - expected1) > 1.0e-20_RKIND) &
         stop "FAIL: horizontal MMS first-level source mismatch"
      if (abs(rho_charge(2,1) - expected2) > 1.0e-20_RKIND) &
         stop "FAIL: horizontal MMS top-level source mismatch"
   end subroutine test_horizontal_mms_source_uses_discrete_vertical_operator

   subroutine test_sphere_mms_source_uses_boundary_compatible_vertical_basis()
      integer, parameter :: nCells = 1, nVertLevels = 1
      real(kind=RKIND), parameter :: epsilon0 = 8.8541878128e-12_RKIND
      real(kind=RKIND), parameter :: sphere_radius = 6371229.0_RKIND
      real(kind=RKIND), parameter :: z_top = 20000.0_RKIND
      real(kind=RKIND) :: latCell(nCells), lonCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND) :: horizontal_eigenvalue, kz, phi_exact, y42, expected

      latCell = [0.25_RKIND]
      lonCell = [0.5_RKIND]
      zMid(:,1) = [0.3_RKIND * z_top]

      call electrostatic_fill_mms_sphere_source(nCells, nVertLevels, latCell, lonCell, zMid, &
                                                epsilon0, sphere_radius, z_top, rho_charge)

      y42 = cos(latCell(1))**2 * (7.0_RKIND * sin(latCell(1))**2 - 1.0_RKIND) * &
            cos(2.0_RKIND * lonCell(1))
      kz = 0.5_RKIND * acos(-1.0_RKIND) / z_top
      horizontal_eigenvalue = 20.0_RKIND / (sphere_radius * sphere_radius)
      phi_exact = y42 * sin(kz * zMid(1,1))
      expected = epsilon0 * (horizontal_eigenvalue + kz*kz) * phi_exact

      if (abs(rho_charge(1,1) - expected) > 1.0e-24_RKIND) &
         stop "FAIL: spherical MMS source mismatch"
   end subroutine test_sphere_mms_source_uses_boundary_compatible_vertical_basis

   subroutine test_sphere_horizontal_mms_source_uses_discrete_vertical_operator()
      integer, parameter :: nCells = 1, nVertLevels = 2
      real(kind=RKIND), parameter :: epsilon0 = 8.8541878128e-12_RKIND
      real(kind=RKIND), parameter :: sphere_radius = 6371229.0_RKIND
      real(kind=RKIND), parameter :: z_top = 20000.0_RKIND
      real(kind=RKIND) :: latCell(nCells), lonCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND) :: v_weight_upper(nVertLevels,nCells)
      real(kind=RKIND) :: v_weight_lower(nVertLevels,nCells)
      real(kind=RKIND) :: volume(nVertLevels,nCells)
      real(kind=RKIND) :: horizontal_coeff, phi1, phi2, y42, kz, expected1, expected2

      latCell = [0.25_RKIND]
      lonCell = [0.5_RKIND]
      zMid(:,1) = [0.25_RKIND * z_top, 0.75_RKIND * z_top]
      volume(:,1) = [2.0_RKIND, 3.0_RKIND]
      v_weight_upper(:,1) = [5.0_RKIND, 0.0_RKIND]
      v_weight_lower(:,1) = [7.0_RKIND, 11.0_RKIND]

      call electrostatic_fill_mms_sphere_horizontal_source(nCells, nVertLevels, latCell, lonCell, &
                                                           zMid, epsilon0, sphere_radius, z_top, &
                                                           v_weight_upper, v_weight_lower, volume, &
                                                           rho_charge)

      y42 = cos(latCell(1))**2 * (7.0_RKIND * sin(latCell(1))**2 - 1.0_RKIND) * &
            cos(2.0_RKIND * lonCell(1))
      kz = 0.5_RKIND * acos(-1.0_RKIND) / z_top
      horizontal_coeff = epsilon0 * 20.0_RKIND / (sphere_radius * sphere_radius)
      phi1 = y42 * sin(kz * zMid(1,1))
      phi2 = y42 * sin(kz * zMid(2,1))

      expected1 = horizontal_coeff * phi1 + &
                  (v_weight_upper(1,1) * (phi1 - phi2) + &
                   v_weight_lower(1,1) * phi1) / volume(1,1)
      expected2 = horizontal_coeff * phi2 + &
                  (v_weight_lower(2,1) * (phi2 - phi1)) / volume(2,1)

      if (abs(rho_charge(1,1) - expected1) > 1.0e-24_RKIND) &
         stop "FAIL: spherical horizontal MMS first-level source mismatch"
      if (abs(rho_charge(2,1) - expected2) > 1.0e-24_RKIND) &
         stop "FAIL: spherical horizontal MMS top-level source mismatch"
   end subroutine test_sphere_horizontal_mms_source_uses_discrete_vertical_operator

   subroutine test_dipole_source_uses_explicit_center_with_expected_lobe_signs()
      integer, parameter :: nCells = 3, nVertLevels = 3
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND) :: center_column_peak, off_center_peak
      real(kind=RKIND), parameter :: x0 = 42000.0_RKIND
      real(kind=RKIND), parameter :: y0 = 42000.0_RKIND
      integer :: iCell

      xCell = [42000.0_RKIND, 84000.0_RKIND, 126000.0_RKIND]
      yCell = [42000.0_RKIND, 42000.0_RKIND, 42000.0_RKIND]
      do iCell = 1, nCells
         zMid(:,iCell) = [5000.0_RKIND, 9000.0_RKIND, 13000.0_RKIND]
      end do

      call electrostatic_fill_dipole_source(nCells, nVertLevels, xCell, yCell, zMid, x0, y0, rho_charge)

      center_column_peak = maxval(abs(rho_charge(:,1)))
      off_center_peak = max(maxval(abs(rho_charge(:,2))), maxval(abs(rho_charge(:,3))))
      if (center_column_peak <= off_center_peak) stop "FAIL: dipole source does not honor explicit center"

      if (rho_charge(1,1) <= 0.0_RKIND) stop "FAIL: lower dipole lobe should be positive"
      if (rho_charge(2,1) >= 0.0_RKIND) stop "FAIL: upper dipole lobe should be negative"
   end subroutine test_dipole_source_uses_explicit_center_with_expected_lobe_signs

   subroutine test_tripole_source_uses_explicit_center_with_expected_lobe_signs()
      integer, parameter :: nCells = 3, nVertLevels = 5
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND) :: center_column_peak, off_center_peak
      real(kind=RKIND), parameter :: x0 = 42000.0_RKIND
      real(kind=RKIND), parameter :: y0 = 42000.0_RKIND
      integer :: iCell

      xCell = [42000.0_RKIND, 84000.0_RKIND, 126000.0_RKIND]
      yCell = [42000.0_RKIND, 42000.0_RKIND, 42000.0_RKIND]
      do iCell = 1, nCells
         zMid(:,iCell) = [4000.0_RKIND, 7000.0_RKIND, 10000.0_RKIND, 13000.0_RKIND, 16000.0_RKIND]
      end do

      call electrostatic_fill_tripole_source(nCells, nVertLevels, xCell, yCell, zMid, x0, y0, rho_charge)

      center_column_peak = maxval(abs(rho_charge(:,1)))
      off_center_peak = max(maxval(abs(rho_charge(:,2))), maxval(abs(rho_charge(:,3))))
      if (center_column_peak <= off_center_peak) stop "FAIL: tripole source does not honor explicit center"

      if (rho_charge(1,1) <= 0.0_RKIND) stop "FAIL: lower tripole lobe should be positive"
      if (rho_charge(2,1) >= 0.0_RKIND) stop "FAIL: middle tripole lobe should be negative"
      if (rho_charge(3,1) <= 0.0_RKIND) stop "FAIL: upper tripole lobe should be positive"
   end subroutine test_tripole_source_uses_explicit_center_with_expected_lobe_signs

   subroutine test_gaussian_monopole_source_is_centered_and_positive()
      integer, parameter :: nCells = 3, nVertLevels = 3
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND), parameter :: x0 = 42000.0_RKIND
      real(kind=RKIND), parameter :: y0 = 42000.0_RKIND
      real(kind=RKIND), parameter :: z0 = 7000.0_RKIND
      integer :: iCell

      xCell = [42000.0_RKIND, 50000.0_RKIND, 58000.0_RKIND]
      yCell = [42000.0_RKIND, 42000.0_RKIND, 42000.0_RKIND]
      do iCell = 1, nCells
         zMid(:,iCell) = [5000.0_RKIND, 7000.0_RKIND, 9000.0_RKIND]
      end do

      call electrostatic_fill_gaussian_monopole_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                       x0, y0, z0, 20.0_RKIND, 2000.0_RKIND, &
                                                       rho_charge)

      if (rho_charge(2,1) <= 0.0_RKIND) stop "FAIL: Gaussian monopole center should be positive"
      if (rho_charge(2,1) <= rho_charge(2,2)) stop "FAIL: Gaussian monopole does not honor explicit center"
      if (any(rho_charge < 0.0_RKIND)) stop "FAIL: Gaussian monopole should not contain negative charge"
   end subroutine test_gaussian_monopole_source_is_centered_and_positive

   subroutine test_gaussian_dipole_sources_have_expected_signs()
      integer, parameter :: nCells = 3, nVertLevels = 5
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_y(nVertLevels,nCells), rho_z(nVertLevels,nCells)
      real(kind=RKIND), parameter :: x0 = 42000.0_RKIND
      real(kind=RKIND), parameter :: y0 = 42000.0_RKIND
      real(kind=RKIND), parameter :: z0 = 8000.0_RKIND
      integer :: iCell

      xCell = [42000.0_RKIND, 42000.0_RKIND, 42000.0_RKIND]
      yCell = [38000.0_RKIND, 42000.0_RKIND, 46000.0_RKIND]
      do iCell = 1, nCells
         zMid(:,iCell) = [4000.0_RKIND, 6000.0_RKIND, 8000.0_RKIND, 10000.0_RKIND, 12000.0_RKIND]
      end do

      call electrostatic_fill_gaussian_dipole_y_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                       x0, y0, z0, 20.0_RKIND, 1200.0_RKIND, &
                                                       8000.0_RKIND, rho_y)
      call electrostatic_fill_gaussian_dipole_z_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                       x0, y0, z0, 20.0_RKIND, 1200.0_RKIND, &
                                                       4000.0_RKIND, rho_z)

      if (rho_y(3,1) <= 0.0_RKIND) stop "FAIL: horizontal dipole positive lobe sign mismatch"
      if (rho_y(3,3) >= 0.0_RKIND) stop "FAIL: horizontal dipole negative lobe sign mismatch"
      if (rho_z(2,2) <= 0.0_RKIND) stop "FAIL: vertical dipole lower positive lobe sign mismatch"
      if (rho_z(4,2) >= 0.0_RKIND) stop "FAIL: vertical dipole upper negative lobe sign mismatch"
   end subroutine test_gaussian_dipole_sources_have_expected_signs

   subroutine test_gaussian_charges_have_expected_integrals()
      integer, parameter :: nSide = 25, nCells = nSide * nSide, nVertLevels = 25
      real(kind=RKIND) :: xCell(nCells), yCell(nCells), zMid(nVertLevels,nCells)
      real(kind=RKIND) :: rho_charge(nVertLevels,nCells)
      real(kind=RKIND), parameter :: x0 = 0.0_RKIND
      real(kind=RKIND), parameter :: y0 = 0.0_RKIND
      real(kind=RKIND), parameter :: z0 = 0.0_RKIND
      real(kind=RKIND), parameter :: charge = 20.0_RKIND
      real(kind=RKIND), parameter :: sigma = 1000.0_RKIND
      real(kind=RKIND), parameter :: separation = 7000.0_RKIND
      real(kind=RKIND), parameter :: half_width = 10000.0_RKIND
      real(kind=RKIND) :: dx, dy, dz, volume, integral
      integer :: i, j, k, iCell

      dx = 2.0_RKIND * half_width / real(nSide, kind=RKIND)
      dy = dx
      dz = 2.0_RKIND * half_width / real(nVertLevels, kind=RKIND)
      volume = dx * dy * dz

      do j = 1, nSide
         do i = 1, nSide
            iCell = (j - 1) * nSide + i
            xCell(iCell) = x0 - half_width + (real(i, kind=RKIND) - 0.5_RKIND) * dx
            yCell(iCell) = y0 - half_width + (real(j, kind=RKIND) - 0.5_RKIND) * dy
         end do
      end do

      do iCell = 1, nCells
         do k = 1, nVertLevels
            zMid(k,iCell) = z0 - half_width + (real(k, kind=RKIND) - 0.5_RKIND) * dz
         end do
      end do

      call electrostatic_fill_gaussian_monopole_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                       x0, y0, z0, charge, sigma, rho_charge)
      integral = sum(rho_charge) * volume
      if (abs(integral - charge) > 1.0e-2_RKIND * charge) &
         stop "FAIL: Gaussian monopole integrated charge mismatch"

      call electrostatic_fill_gaussian_dipole_y_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                       x0, y0, z0, charge, sigma, separation, rho_charge)
      integral = sum(rho_charge) * volume
      if (abs(integral) > 1.0e-10_RKIND * charge) &
         stop "FAIL: horizontal Gaussian dipole should be net neutral"

      call electrostatic_fill_gaussian_dipole_z_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                       x0, y0, z0, charge, sigma, separation, rho_charge)
      integral = sum(rho_charge) * volume
      if (abs(integral) > 1.0e-10_RKIND * charge) &
         stop "FAIL: vertical Gaussian dipole should be net neutral"
   end subroutine test_gaussian_charges_have_expected_integrals

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
