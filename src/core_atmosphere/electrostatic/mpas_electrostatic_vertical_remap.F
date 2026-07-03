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
   public :: electrostatic_build_poisson_grid
   public :: electrostatic_build_poisson_fv_weights

contains

!***********************************************************************
!
!  routine electrostatic_build_poisson_grid
!
!> \brief Build the uniform Poisson z grid and per-column active air bands
!
!-----------------------------------------------------------------------
   subroutine electrostatic_build_poisson_grid(nCells, nVertLevels, zGround, z_bottom, z_top, &
                                               zFace, poissonZmid, active, kFirstActive, airThickness)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: zGround
      real(kind=RKIND), intent(in) :: z_bottom
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nVertLevels+1), intent(out) :: zFace
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: poissonZmid
      logical, dimension(nVertLevels,nCells), intent(out) :: active
      integer, dimension(nCells), intent(out) :: kFirstActive
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: airThickness

      integer :: iCell, k, kf
      real(kind=RKIND) :: dzP, t
      real(kind=RKIND), parameter :: sliver_frac = 0.1_RKIND

      dzP = (z_top - z_bottom) / real(nVertLevels, kind=RKIND)
      do k = 1, nVertLevels + 1
         zFace(k) = z_bottom + real(k-1, kind=RKIND) * dzP
      end do

      do iCell = 1, nCells
         kf = nVertLevels + 1
         do k = 1, nVertLevels
            t = max(0.0_RKIND, zFace(k+1) - max(zFace(k), zGround(iCell)))
            if (t > 0.0_RKIND) then
               kf = k
               exit
            end if
         end do

         if (kf <= nVertLevels) then
            if (zFace(kf+1) - zGround(iCell) < sliver_frac * dzP .and. kf < nVertLevels) kf = kf + 1
         end if
         kFirstActive(iCell) = kf

         do k = 1, nVertLevels
            if (k < kf) then
               active(k,iCell) = .false.
               airThickness(k,iCell) = 0.0_RKIND
               poissonZmid(k,iCell) = 0.0_RKIND
            else if (k == kf) then
               active(k,iCell) = .true.
               airThickness(k,iCell) = zFace(k+1) - zGround(iCell)
               poissonZmid(k,iCell) = 0.5_RKIND * (zGround(iCell) + zFace(k+1))
            else
               active(k,iCell) = .true.
               airThickness(k,iCell) = zFace(k+1) - zFace(k)
               poissonZmid(k,iCell) = 0.5_RKIND * (zFace(k) + zFace(k+1))
            end if
         end do
      end do

   end subroutine electrostatic_build_poisson_grid

!***********************************************************************
!
!  routine electrostatic_build_poisson_fv_weights
!
!> \brief Build terrain-aware finite-volume operator weights
!
!-----------------------------------------------------------------------
   subroutine electrostatic_build_poisson_fv_weights(nCells, nEdges, nVertLevels, cellsOnEdge, &
                                                     areaCell, dcEdge, dvEdge, epsilon0, &
                                                     zFace, zGround, z_top, poissonZmid, &
                                                     active, kFirstActive, airThickness, &
                                                     h_weight, v_weight_upper, v_weight_lower, &
                                                     diag_extra, volume)

      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
      real(kind=RKIND), dimension(nCells), intent(in) :: areaCell
      real(kind=RKIND), dimension(nEdges), intent(in) :: dcEdge
      real(kind=RKIND), dimension(nEdges), intent(in) :: dvEdge
      real(kind=RKIND), intent(in) :: epsilon0
      real(kind=RKIND), dimension(nVertLevels+1), intent(in) :: zFace
      real(kind=RKIND), dimension(nCells), intent(in) :: zGround
      real(kind=RKIND), intent(in) :: z_top
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: poissonZmid
      logical, dimension(nVertLevels,nCells), intent(in) :: active
      integer, dimension(nCells), intent(in) :: kFirstActive
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: airThickness
      real(kind=RKIND), dimension(nVertLevels,nEdges), intent(out) :: h_weight
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: v_weight_upper
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: v_weight_lower
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: diag_extra
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: volume

      integer :: iCell, iEdge, k, c1, c2
      real(kind=RKIND) :: lo, hi, tOpen, tW1, tW2, wall_scale

      h_weight = 0.0_RKIND
      v_weight_upper = 0.0_RKIND
      v_weight_lower = 0.0_RKIND
      diag_extra = 0.0_RKIND
      volume = 0.0_RKIND

      do iCell = 1, nCells
         do k = 1, nVertLevels
            volume(k,iCell) = areaCell(iCell) * airThickness(k,iCell)
         end do

         do k = kFirstActive(iCell), nVertLevels - 1
            v_weight_upper(k,iCell) = epsilon0 * areaCell(iCell) / &
                                      (poissonZmid(k+1,iCell) - poissonZmid(k,iCell))
            v_weight_lower(k+1,iCell) = v_weight_upper(k,iCell)
         end do

         if (kFirstActive(iCell) <= nVertLevels) then
            k = kFirstActive(iCell)
            v_weight_lower(k,iCell) = epsilon0 * areaCell(iCell) / &
                                      (poissonZmid(k,iCell) - zGround(iCell))
         end if
      end do

      do iEdge = 1, nEdges
         c1 = cellsOnEdge(1,iEdge)
         c2 = cellsOnEdge(2,iEdge)
         if (c1 < 1 .or. c1 > nCells .or. c2 < 1 .or. c2 > nCells) cycle
         wall_scale = epsilon0 * dvEdge(iEdge) / (0.5_RKIND * dcEdge(iEdge))
         do k = 1, nVertLevels
            lo = zFace(k)
            hi = min(zFace(k+1), z_top)
            if (hi <= lo) cycle
            tOpen = max(0.0_RKIND, hi - max(lo, zGround(c1), zGround(c2)))
            tW1 = max(0.0_RKIND, min(hi, zGround(c2)) - max(lo, zGround(c1)))
            tW2 = max(0.0_RKIND, min(hi, zGround(c1)) - max(lo, zGround(c2)))
            if (active(k,c1) .and. active(k,c2)) then
               h_weight(k,iEdge) = epsilon0 * dvEdge(iEdge) * tOpen / dcEdge(iEdge)
            else if (active(k,c1)) then
               tW1 = tW1 + tOpen
            else if (active(k,c2)) then
               tW2 = tW2 + tOpen
            end if
            if (active(k,c1)) diag_extra(k,c1) = diag_extra(k,c1) + wall_scale * tW1
            if (active(k,c2)) diag_extra(k,c2) = diag_extra(k,c2) + wall_scale * tW2
         end do
      end do

   end subroutine electrostatic_build_poisson_fv_weights

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
