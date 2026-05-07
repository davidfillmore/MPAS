! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_electrostatic_benchmarks
!
!> \brief Analytic benchmark utilities for MPAS electrostatics
!
!-----------------------------------------------------------------------
module mpas_electrostatic_benchmarks

   use mpas_kind_types, only : RKIND

   implicit none
   private
   save

   public :: electrostatic_mms_cart_phi_exact
   public :: electrostatic_compute_error_norms

contains

!***********************************************************************
!
!  routine electrostatic_mms_cart_phi_exact
!
!> \brief Return the corrected Cartesian MMS electrostatic potential
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function electrostatic_mms_cart_phi_exact(x, y, z, x_period, y_period, z_top)

      real(kind=RKIND), intent(in) :: x
      real(kind=RKIND), intent(in) :: y
      real(kind=RKIND), intent(in) :: z
      real(kind=RKIND), intent(in) :: x_period
      real(kind=RKIND), intent(in) :: y_period
      real(kind=RKIND), intent(in) :: z_top

      real(kind=RKIND) :: kx, ky, kz, pi

      pi = acos(-1.0_RKIND)
      kx = 2.0_RKIND * pi / x_period
      ky = 2.0_RKIND * pi / y_period
      kz = 0.5_RKIND * pi / z_top

      electrostatic_mms_cart_phi_exact = sin(kx*x) * sin(ky*y) * sin(kz*z)

   end function electrostatic_mms_cart_phi_exact

!***********************************************************************
!
!  routine electrostatic_compute_error_norms
!
!> \brief Compute relative L2 and absolute Linf error norms for potential
!
!-----------------------------------------------------------------------
   subroutine electrostatic_compute_error_norms(phi, phi_exact, volume, l2, linf)

      real(kind=RKIND), dimension(:,:), intent(in) :: phi
      real(kind=RKIND), dimension(:,:), intent(in) :: phi_exact
      real(kind=RKIND), dimension(:,:), intent(in) :: volume
      real(kind=RKIND), intent(out) :: l2
      real(kind=RKIND), intent(out) :: linf

      real(kind=RKIND) :: denominator, numerator

      numerator = sum(volume * (phi - phi_exact)**2)
      denominator = sum(volume * phi_exact**2)

      if (denominator > 0.0_RKIND) then
         l2 = sqrt(numerator / denominator)
      else if (numerator > 0.0_RKIND) then
         l2 = sqrt(numerator)
      else
         l2 = 0.0_RKIND
      end if

      linf = maxval(abs(phi - phi_exact))

   end subroutine electrostatic_compute_error_norms

end module mpas_electrostatic_benchmarks
