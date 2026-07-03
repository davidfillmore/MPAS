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
   public :: electrostatic_horizontal_metric_scale

contains

!***********************************************************************
!
!  routine electrostatic_apply_ground_rhs
!
!> \brief Add ground Dirichlet potential contributions to the RHS
!
!-----------------------------------------------------------------------
   subroutine electrostatic_apply_ground_rhs(nCells, nVertLevels, phi_ground, v_weight_lower, rhs, &
                                             kFirstActive, wall_weight)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), intent(in) :: phi_ground
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_lower
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: rhs
      integer, dimension(nCells), intent(in), optional :: kFirstActive
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in), optional :: wall_weight

      integer :: iCell, k, kf

      do iCell = 1, nCells
         kf = 1
         if (present(kFirstActive)) kf = kFirstActive(iCell)
         if (kf > nVertLevels) cycle   ! fully-buried column sentinel (kf = nVertLevels+1)

         rhs(kf,iCell) = rhs(kf,iCell) + v_weight_lower(kf,iCell) * phi_ground

         if (present(wall_weight)) then
            ! Grounded lateral walls: the wall Dirichlet sink on the diagonal
            ! (diag_extra) must be paired with the same datum on the RHS, or
            ! the walls are pinned to 0 while the bottom is pinned to
            ! phi_ground. Builder contract: wall_weight(k,i) == 0 on inactive
            ! levels, so no separate active mask is needed here.
            do k = 1, nVertLevels
               rhs(k,iCell) = rhs(k,iCell) + wall_weight(k,iCell) * phi_ground
            end do
         end if
      end do

   end subroutine electrostatic_apply_ground_rhs

!***********************************************************************
!
!  routine electrostatic_compute_vertical_E
!
!> \brief Compute cell-centered vertical electric field from potential
!
!-----------------------------------------------------------------------
   subroutine electrostatic_compute_vertical_E(nCells, nVertLevels, phi, zMid, phi_ground, z_top, &
                                               zGround, kFirstActive, E_vertical)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: phi
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: zMid
      real(kind=RKIND), intent(in) :: phi_ground
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nCells), intent(in) :: zGround
      integer, dimension(nCells), intent(in) :: kFirstActive
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: E_vertical

      integer :: iCell, k, kf
      real(kind=RKIND) :: z0, z1, z2, denom

      E_vertical = 0.0_RKIND

      do iCell = 1, nCells
         kf = kFirstActive(iCell)
         if (kf > nVertLevels) cycle

         if (kf >= nVertLevels) then
            E_vertical(kf,iCell) = -(phi(kf,iCell) - phi_ground) / &
                                   (zMid(kf,iCell) - zGround(iCell))
            cycle
         end if

         z0 = zGround(iCell)
         z1 = zMid(kf,iCell)
         z2 = zMid(kf+1,iCell)
         E_vertical(kf,iCell) = -( &
            phi_ground      * (z1 - z2) / ((z0 - z1) * (z0 - z2)) + &
            phi(kf,iCell)   * (2.0_RKIND * z1 - z0 - z2) / ((z1 - z0) * (z1 - z2)) + &
            phi(kf+1,iCell) * (z1 - z0) / ((z2 - z0) * (z2 - z1)) )

         do k = kf + 1, nVertLevels - 1
            E_vertical(k,iCell) = -(phi(k+1,iCell) - phi(k-1,iCell)) / &
                                   (zMid(k+1,iCell) - zMid(k-1,iCell))
         end do

         z1 = zMid(nVertLevels-1,iCell)
         z2 = zMid(nVertLevels,iCell)
         denom = (z1 - z_top)**2 - (z2 - z_top)**2
         E_vertical(nVertLevels,iCell) = &
            -2.0_RKIND * (phi(nVertLevels-1,iCell) - phi(nVertLevels,iCell)) * &
            (z2 - z_top) / denom
      end do

   end subroutine electrostatic_compute_vertical_E

!***********************************************************************
!
!  routine electrostatic_horizontal_metric_scale
!
!> \brief Return the scale factor needed for horizontal mesh metrics
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function electrostatic_horizontal_metric_scale(nEdges, dcEdge, sphere_radius, on_a_sphere)

      integer, intent(in) :: nEdges
      real(kind=RKIND), dimension(nEdges), intent(in) :: dcEdge
      real(kind=RKIND), intent(in) :: sphere_radius
      logical, intent(in) :: on_a_sphere

      ! Detect whether the horizontal mesh metrics are dimensional (metres,
      ! the standard init_atmosphere output) or non-dimensional (radians, a
      ! raw unit-sphere SCVT mesh). The separation is unambiguous: no physical
      ! atmospheric mesh has cell spacing below ~1 m, and no unit-sphere
      ! angular spacing exceeds pi rad, so a maximum |dcEdge| below this floor
      ! marks a non-dimensional mesh whose metrics must be scaled to metres by
      ! the configured sphere radius (arc length = R * angle).
      real(kind=RKIND), parameter :: dimensional_metric_floor = 10.0_RKIND  ! metres

      electrostatic_horizontal_metric_scale = 1.0_RKIND

      if (.not. on_a_sphere) return
      if (nEdges <= 0) return
      if (sphere_radius <= 1.0_RKIND) return

      if (maxval(abs(dcEdge(1:nEdges))) < dimensional_metric_floor) then
         electrostatic_horizontal_metric_scale = sphere_radius
      end if

   end function electrostatic_horizontal_metric_scale

end module mpas_electrostatic_bcs
