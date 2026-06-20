! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_electrostatic_vertical_remap
!
!> \brief Vertical interpolation utilities for MPAS electrostatic fields
!
!-----------------------------------------------------------------------
module mpas_electrostatic_vertical_remap

   use mpas_kind_types, only : RKIND

   implicit none
   private
   save

   public :: electrostatic_interp_mpas_to_poisson
   public :: electrostatic_interp_w_to_poisson_mid

contains

!***********************************************************************
!
!  routine electrostatic_interp_mpas_to_poisson
!
!> \brief Interpolate cell-centered MPAS fields to Poisson altitude levels
!
!-----------------------------------------------------------------------
   subroutine electrostatic_interp_mpas_to_poisson(nCells, nVertLevels, &
                                                   zMpasMid, zPoissonMid, active, source, target)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMpasMid
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zPoissonMid
      logical, dimension(nVertLevels,nCells), intent(in) :: active
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: source
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: target

      integer :: iCell, k

      target = 0.0_RKIND

      do iCell = 1, nCells
         do k = 1, nVertLevels
            if (active(k,iCell)) then
               target(k,iCell) = interp_column_linear(nVertLevels, zMpasMid(:,iCell), &
                                                       source(:,iCell), zPoissonMid(k,iCell))
            end if
         end do
      end do

   end subroutine electrostatic_interp_mpas_to_poisson

!***********************************************************************
!
!  routine electrostatic_interp_w_to_poisson_mid
!
!> \brief Interpolate MPAS w-interface values to Poisson midpoints
!
!-----------------------------------------------------------------------
   subroutine electrostatic_interp_w_to_poisson_mid(nCells, nVertLevels, &
                                                    zgrid, zPoissonMid, active, w, wPoissonMid)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels+1,nCells), intent(in) :: zgrid
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zPoissonMid
      logical, dimension(nVertLevels,nCells), intent(in) :: active
      real(kind=RKIND), dimension(nVertLevels+1,nCells), intent(in) :: w
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: wPoissonMid

      integer :: iCell, k

      wPoissonMid = 0.0_RKIND

      do iCell = 1, nCells
         do k = 1, nVertLevels
            if (active(k,iCell)) then
               wPoissonMid(k,iCell) = interp_column_linear(nVertLevels+1, zgrid(:,iCell), &
                                                            w(:,iCell), zPoissonMid(k,iCell))
            end if
         end do
      end do

   end subroutine electrostatic_interp_w_to_poisson_mid

!***********************************************************************
!
!  function interp_column_linear
!
!> \brief Piecewise-linear interpolation with endpoint clamping
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function interp_column_linear(nSrc, zSrc, fSrc, zTarget)

      integer, intent(in) :: nSrc
      real(kind=RKIND), dimension(nSrc), intent(in) :: zSrc
      real(kind=RKIND), dimension(nSrc), intent(in) :: fSrc
      real(kind=RKIND), intent(in) :: zTarget

      integer :: k
      real(kind=RKIND) :: weight

      if (nSrc <= 1) then
         interp_column_linear = fSrc(1)
         return
      end if

      if (zTarget <= zSrc(1)) then
         interp_column_linear = fSrc(1)
         return
      end if

      if (zTarget >= zSrc(nSrc)) then
         interp_column_linear = fSrc(nSrc)
         return
      end if

      do k = 1, nSrc - 1
         if (zTarget >= zSrc(k) .and. zTarget <= zSrc(k+1)) then
            if (zSrc(k+1) <= zSrc(k)) then
               interp_column_linear = fSrc(k)
            else
               weight = (zTarget - zSrc(k)) / (zSrc(k+1) - zSrc(k))
               interp_column_linear = (1.0_RKIND - weight) * fSrc(k) + weight * fSrc(k+1)
            end if
            return
         end if
      end do

      interp_column_linear = fSrc(nSrc)

   end function interp_column_linear

end module mpas_electrostatic_vertical_remap
