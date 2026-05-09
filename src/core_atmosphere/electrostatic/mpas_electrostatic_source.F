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
