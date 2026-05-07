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

end module mpas_electrostatic_source
