! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_electrostatic_source
!
!> \brief Charge-density source utilities for MPAS electrostatics
!
!-----------------------------------------------------------------------
module mpas_electrostatic_source

   use mpas_kind_types, only : RKIND

   implicit none
   private
   save

   public :: electrostatic_fill_mms_cart_source
   public :: electrostatic_fill_mms_cart_horizontal_source
   public :: electrostatic_fill_mms_sphere_source
   public :: electrostatic_fill_mms_sphere_horizontal_source
   public :: electrostatic_fill_dipole_source
   public :: electrostatic_fill_tripole_source
   public :: electrostatic_fill_gaussian_monopole_source
   public :: electrostatic_fill_gaussian_dipole_y_source
   public :: electrostatic_fill_gaussian_dipole_z_source

contains

!***********************************************************************
!
!  routine electrostatic_fill_mms_cart_source
!
!> \brief Fill charge density for the corrected Cartesian MMS case
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_mms_cart_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                 epsilon0, x_period, y_period, z_top, rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: epsilon0
      real(kind=RKIND), intent(in) :: x_period
      real(kind=RKIND), intent(in) :: y_period
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      integer :: iCell, k
      real(kind=RKIND) :: coeff, kx, ky, kz, phi_exact, pi

      pi = acos(-1.0_RKIND)
      kx = 2.0_RKIND * pi / x_period
      ky = 2.0_RKIND * pi / y_period
      kz = 0.5_RKIND * pi / z_top
      coeff = epsilon0 * (kx*kx + ky*ky + kz*kz)

      do iCell = 1, nCells
         do k = 1, nVertLevels
            phi_exact = sin(kx*xCell(iCell)) * sin(ky*yCell(iCell)) * sin(kz*zMid(k,iCell))
            rho_charge(k,iCell) = coeff * phi_exact
         end do
      end do

   end subroutine electrostatic_fill_mms_cart_source

!***********************************************************************
!
!  routine electrostatic_fill_mms_cart_horizontal_source
!
!> \brief Fill a Cartesian MMS source that isolates horizontal truncation error
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_mms_cart_horizontal_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                            epsilon0, x_period, y_period, z_top, &
                                                            v_weight_upper, v_weight_lower, volume, &
                                                            rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: epsilon0
      real(kind=RKIND), intent(in) :: x_period
      real(kind=RKIND), intent(in) :: y_period
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_upper
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_lower
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: volume
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      integer :: iCell, k
      real(kind=RKIND) :: horizontal_coeff, kx, ky, kz, phi_here, phi_lower, phi_upper, pi
      real(kind=RKIND) :: vertical_rhs

      pi = acos(-1.0_RKIND)
      kx = 2.0_RKIND * pi / x_period
      ky = 2.0_RKIND * pi / y_period
      kz = 0.5_RKIND * pi / z_top
      horizontal_coeff = epsilon0 * (kx*kx + ky*ky)

      do iCell = 1, nCells
         do k = 1, nVertLevels
            phi_here = sin(kx*xCell(iCell)) * sin(ky*yCell(iCell)) * sin(kz*zMid(k,iCell))

            if (k < nVertLevels) then
               phi_upper = sin(kx*xCell(iCell)) * sin(ky*yCell(iCell)) * sin(kz*zMid(k+1,iCell))
            else
               phi_upper = phi_here
            end if

            if (k > 1) then
               phi_lower = sin(kx*xCell(iCell)) * sin(ky*yCell(iCell)) * sin(kz*zMid(k-1,iCell))
            else
               phi_lower = 0.0_RKIND
            end if

            vertical_rhs = 0.0_RKIND
            if (k < nVertLevels) then
               vertical_rhs = vertical_rhs + v_weight_upper(k,iCell) * (phi_here - phi_upper)
            end if
            vertical_rhs = vertical_rhs + v_weight_lower(k,iCell) * (phi_here - phi_lower)

            rho_charge(k,iCell) = horizontal_coeff * phi_here + vertical_rhs / volume(k,iCell)
         end do
      end do

   end subroutine electrostatic_fill_mms_cart_horizontal_source

!***********************************************************************
!
!  routine electrostatic_fill_mms_sphere_source
!
!> \brief Fill charge density for the spherical-harmonic MMS case
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_mms_sphere_source(nCells, nVertLevels, latCell, lonCell, zMid, &
                                                   epsilon0, sphere_radius, z_top, rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: latCell
      real(kind=RKIND), dimension(nCells), intent(in) :: lonCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: epsilon0
      real(kind=RKIND), intent(in) :: sphere_radius
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      integer :: iCell, k
      real(kind=RKIND) :: coeff, horizontal_eigenvalue, kz, radius

      radius = max(sphere_radius, 1.0_RKIND)
      kz = 0.5_RKIND * acos(-1.0_RKIND) / z_top
      horizontal_eigenvalue = 20.0_RKIND / (radius * radius)
      coeff = epsilon0 * (horizontal_eigenvalue + kz*kz)

      do iCell = 1, nCells
         do k = 1, nVertLevels
            rho_charge(k,iCell) = coeff * sphere_mms_phi(latCell(iCell), lonCell(iCell), zMid(k,iCell), z_top)
         end do
      end do

   end subroutine electrostatic_fill_mms_sphere_source

!***********************************************************************
!
!  routine electrostatic_fill_mms_sphere_horizontal_source
!
!> \brief Fill a spherical MMS source that isolates horizontal truncation error
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_mms_sphere_horizontal_source(nCells, nVertLevels, latCell, lonCell, zMid, &
                                                              epsilon0, sphere_radius, z_top, &
                                                              v_weight_upper, v_weight_lower, volume, &
                                                              rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: latCell
      real(kind=RKIND), dimension(nCells), intent(in) :: lonCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: epsilon0
      real(kind=RKIND), intent(in) :: sphere_radius
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_upper
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_lower
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: volume
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      integer :: iCell, k
      real(kind=RKIND) :: horizontal_coeff, phi_here, phi_lower, phi_upper, radius, vertical_rhs

      radius = max(sphere_radius, 1.0_RKIND)
      horizontal_coeff = epsilon0 * 20.0_RKIND / (radius * radius)

      do iCell = 1, nCells
         do k = 1, nVertLevels
            phi_here = sphere_mms_phi(latCell(iCell), lonCell(iCell), zMid(k,iCell), z_top)

            if (k < nVertLevels) then
               phi_upper = sphere_mms_phi(latCell(iCell), lonCell(iCell), zMid(k+1,iCell), z_top)
            else
               phi_upper = phi_here
            end if

            if (k > 1) then
               phi_lower = sphere_mms_phi(latCell(iCell), lonCell(iCell), zMid(k-1,iCell), z_top)
            else
               phi_lower = 0.0_RKIND
            end if

            vertical_rhs = 0.0_RKIND
            if (k < nVertLevels) then
               vertical_rhs = vertical_rhs + v_weight_upper(k,iCell) * (phi_here - phi_upper)
            end if
            vertical_rhs = vertical_rhs + v_weight_lower(k,iCell) * (phi_here - phi_lower)

            rho_charge(k,iCell) = horizontal_coeff * phi_here + vertical_rhs / volume(k,iCell)
         end do
      end do

   end subroutine electrostatic_fill_mms_sphere_horizontal_source

!***********************************************************************
!
!  routine electrostatic_fill_dipole_source
!
!> \brief Fill a neutral synthetic thundercloud dipole charge density
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_dipole_source(nCells, nVertLevels, xCell, yCell, zMid, x0, y0, rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: x0
      real(kind=RKIND), intent(in) :: y0
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      integer :: iCell, k, lobe
      real(kind=RKIND), parameter :: sigma_h = 8000.0_RKIND
      real(kind=RKIND), parameter :: sigma_v = 1500.0_RKIND
      real(kind=RKIND), dimension(2), parameter :: q_lobe = [20.0_RKIND, -20.0_RKIND]
      real(kind=RKIND), dimension(2), parameter :: z_lobe = [5000.0_RKIND, 9000.0_RKIND]
      real(kind=RKIND) :: amp, dx2, dy2, dz2, horizontal_factor, norm, pi

      pi = acos(-1.0_RKIND)
      norm = (2.0_RKIND * pi)**1.5_RKIND * sigma_h * sigma_h * sigma_v

      rho_charge = 0.0_RKIND
      do iCell = 1, nCells
         dx2 = (xCell(iCell) - x0)**2
         dy2 = (yCell(iCell) - y0)**2
         horizontal_factor = exp(-(dx2 + dy2) / (2.0_RKIND * sigma_h * sigma_h))
         do k = 1, nVertLevels
            do lobe = 1, 2
               dz2 = (zMid(k,iCell) - z_lobe(lobe))**2
               amp = q_lobe(lobe) / norm
               rho_charge(k,iCell) = rho_charge(k,iCell) + amp * horizontal_factor * &
                                      exp(-dz2 / (2.0_RKIND * sigma_v * sigma_v))
            end do
         end do
      end do

   end subroutine electrostatic_fill_dipole_source

!***********************************************************************
!
!  routine electrostatic_fill_tripole_source
!
!> \brief Fill a neutral synthetic thundercloud tripole charge density
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_tripole_source(nCells, nVertLevels, xCell, yCell, zMid, x0, y0, rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: x0
      real(kind=RKIND), intent(in) :: y0
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      integer :: iCell, k, lobe
      real(kind=RKIND), parameter :: sigma_h = 8000.0_RKIND
      real(kind=RKIND), parameter :: sigma_v = 1500.0_RKIND
      real(kind=RKIND), dimension(3), parameter :: q_lobe = [20.0_RKIND, -40.0_RKIND, 20.0_RKIND]
      real(kind=RKIND), dimension(3), parameter :: z_lobe = [4000.0_RKIND, 7000.0_RKIND, 10000.0_RKIND]
      real(kind=RKIND) :: amp, dx2, dy2, dz2, horizontal_factor, norm, pi

      pi = acos(-1.0_RKIND)
      norm = (2.0_RKIND * pi)**1.5_RKIND * sigma_h * sigma_h * sigma_v

      rho_charge = 0.0_RKIND
      do iCell = 1, nCells
         dx2 = (xCell(iCell) - x0)**2
         dy2 = (yCell(iCell) - y0)**2
         horizontal_factor = exp(-(dx2 + dy2) / (2.0_RKIND * sigma_h * sigma_h))
         do k = 1, nVertLevels
            do lobe = 1, 3
               dz2 = (zMid(k,iCell) - z_lobe(lobe))**2
               amp = q_lobe(lobe) / norm
               rho_charge(k,iCell) = rho_charge(k,iCell) + amp * horizontal_factor * &
                                      exp(-dz2 / (2.0_RKIND * sigma_v * sigma_v))
            end do
         end do
      end do

   end subroutine electrostatic_fill_tripole_source

!***********************************************************************
!
!  routine electrostatic_fill_gaussian_monopole_source
!
!> \brief Fill an isotropic Gaussian monopole charge density
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_gaussian_monopole_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                          x0, y0, z0, charge, sigma, rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: x0
      real(kind=RKIND), intent(in) :: y0
      real(kind=RKIND), intent(in) :: z0
      real(kind=RKIND), intent(in) :: charge
      real(kind=RKIND), intent(in) :: sigma
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      rho_charge = 0.0_RKIND
      call add_gaussian_charge_lobe(nCells, nVertLevels, xCell, yCell, zMid, x0, y0, z0, &
                                    charge, sigma, rho_charge)

   end subroutine electrostatic_fill_gaussian_monopole_source

!***********************************************************************
!
!  routine electrostatic_fill_gaussian_dipole_y_source
!
!> \brief Fill a y-oriented neutral isotropic Gaussian dipole charge density
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_gaussian_dipole_y_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                          x0, y0, z0, charge, sigma, separation, &
                                                          rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: x0
      real(kind=RKIND), intent(in) :: y0
      real(kind=RKIND), intent(in) :: z0
      real(kind=RKIND), intent(in) :: charge
      real(kind=RKIND), intent(in) :: sigma
      real(kind=RKIND), intent(in) :: separation
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      rho_charge = 0.0_RKIND
      call add_gaussian_charge_lobe(nCells, nVertLevels, xCell, yCell, zMid, &
                                    x0, y0 - 0.5_RKIND * separation, z0, &
                                    charge, sigma, rho_charge)
      call add_gaussian_charge_lobe(nCells, nVertLevels, xCell, yCell, zMid, &
                                    x0, y0 + 0.5_RKIND * separation, z0, &
                                    -charge, sigma, rho_charge)

   end subroutine electrostatic_fill_gaussian_dipole_y_source

!***********************************************************************
!
!  routine electrostatic_fill_gaussian_dipole_z_source
!
!> \brief Fill a z-oriented neutral isotropic Gaussian dipole charge density
!
!-----------------------------------------------------------------------
   subroutine electrostatic_fill_gaussian_dipole_z_source(nCells, nVertLevels, xCell, yCell, zMid, &
                                                          x0, y0, z0, charge, sigma, separation, &
                                                          rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: x0
      real(kind=RKIND), intent(in) :: y0
      real(kind=RKIND), intent(in) :: z0
      real(kind=RKIND), intent(in) :: charge
      real(kind=RKIND), intent(in) :: sigma
      real(kind=RKIND), intent(in) :: separation
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: rho_charge

      rho_charge = 0.0_RKIND
      call add_gaussian_charge_lobe(nCells, nVertLevels, xCell, yCell, zMid, &
                                    x0, y0, z0 - 0.5_RKIND * separation, &
                                    charge, sigma, rho_charge)
      call add_gaussian_charge_lobe(nCells, nVertLevels, xCell, yCell, zMid, &
                                    x0, y0, z0 + 0.5_RKIND * separation, &
                                    -charge, sigma, rho_charge)

   end subroutine electrostatic_fill_gaussian_dipole_z_source

!***********************************************************************
!
!  routine add_gaussian_charge_lobe
!
!> \brief Add one normalized isotropic Gaussian charge lobe
!
!-----------------------------------------------------------------------
   subroutine add_gaussian_charge_lobe(nCells, nVertLevels, xCell, yCell, zMid, &
                                       x0, y0, z0, charge, sigma, rho_charge)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: x0
      real(kind=RKIND), intent(in) :: y0
      real(kind=RKIND), intent(in) :: z0
      real(kind=RKIND), intent(in) :: charge
      real(kind=RKIND), intent(in) :: sigma
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: rho_charge

      integer :: iCell, k
      real(kind=RKIND) :: amp, dx2, dy2, dz2, norm, pi, sigma_safe

      pi = acos(-1.0_RKIND)
      sigma_safe = max(sigma, 1.0_RKIND)
      norm = (2.0_RKIND * pi)**1.5_RKIND * sigma_safe**3
      amp = charge / norm

      do iCell = 1, nCells
         dx2 = (xCell(iCell) - x0)**2
         dy2 = (yCell(iCell) - y0)**2
         do k = 1, nVertLevels
            dz2 = (zMid(k,iCell) - z0)**2
            rho_charge(k,iCell) = rho_charge(k,iCell) + amp * &
                                   exp(-(dx2 + dy2 + dz2) / (2.0_RKIND * sigma_safe**2))
         end do
      end do

   end subroutine add_gaussian_charge_lobe

!***********************************************************************
!
!  routine sphere_mms_phi
!
!> \brief Return the Tier A.2 spherical-harmonic MMS potential
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function sphere_mms_phi(lat, lon, z, z_top)

      real(kind=RKIND), intent(in) :: lat
      real(kind=RKIND), intent(in) :: lon
      real(kind=RKIND), intent(in) :: z
      real(kind=RKIND), intent(in) :: z_top

      real(kind=RKIND) :: kz, y42

      kz = 0.5_RKIND * acos(-1.0_RKIND) / z_top
      y42 = cos(lat)**2 * (7.0_RKIND * sin(lat)**2 - 1.0_RKIND) * cos(2.0_RKIND * lon)
      sphere_mms_phi = y42 * sin(kz * z)

   end function sphere_mms_phi

end module mpas_electrostatic_source
