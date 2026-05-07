! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_electrostatic_bcs
!
!> \brief Boundary-condition helpers for MPAS electrostatics
!
!-----------------------------------------------------------------------
module mpas_electrostatic_bcs

   use mpas_kind_types, only : RKIND

   implicit none
   private
   save

   public :: electrostatic_apply_ground_rhs
   public :: electrostatic_compute_vertical_E

contains

!***********************************************************************
!
!  routine electrostatic_apply_ground_rhs
!
!> \brief Add ground Dirichlet potential contributions to the RHS
!
!-----------------------------------------------------------------------
   subroutine electrostatic_apply_ground_rhs(nCells, nVertLevels, phi_ground, v_weight_lower, rhs)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), intent(in) :: phi_ground
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_lower
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: rhs

      integer :: iCell

      do iCell = 1, nCells
         rhs(1,iCell) = rhs(1,iCell) + v_weight_lower(1,iCell) * phi_ground
      end do

   end subroutine electrostatic_apply_ground_rhs

!***********************************************************************
!
!  routine electrostatic_compute_vertical_E
!
!> \brief Compute cell-centered vertical electric field from potential
!
!-----------------------------------------------------------------------
   subroutine electrostatic_compute_vertical_E(nCells, nVertLevels, phi, zMid, phi_ground, E_vertical)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: phi
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: phi_ground
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: E_vertical

      integer :: iCell, k

      E_vertical = 0.0_RKIND

      do iCell = 1, nCells
         if (nVertLevels == 1) then
            E_vertical(1,iCell) = 0.0_RKIND
         else
            E_vertical(1,iCell) = -(phi(2,iCell) - phi_ground) / zMid(2,iCell)

            do k = 2, nVertLevels - 1
               E_vertical(k,iCell) = -(phi(k+1,iCell) - phi(k-1,iCell)) / &
                                      (zMid(k+1,iCell) - zMid(k-1,iCell))
            end do

            E_vertical(nVertLevels,iCell) = 0.0_RKIND
         end if
      end do

   end subroutine electrostatic_compute_vertical_E

end module mpas_electrostatic_bcs
