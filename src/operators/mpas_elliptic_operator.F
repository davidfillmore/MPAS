! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_elliptic_operator
!
!> \brief MPAS electrostatic elliptic operator kernels
!> \details
!>  Applies the finite-volume operator A = -div(epsilon grad(.)) with a
!>  grounded lower boundary and homogeneous Neumann upper/lateral boundaries.
!
!-----------------------------------------------------------------------
module mpas_elliptic_operator

   use mpas_kind_types, only : RKIND

   implicit none
   private
   save

   public :: mpas_elliptic_operator_type, elliptic_operator_init
   public :: elliptic_operator_init_from_weights

   type :: mpas_elliptic_operator_type
      integer :: nCells = 0
      integer :: nEdges = 0
      integer :: nVertLevels = 0
      integer :: maxEdges = 0
      integer, allocatable :: nEdgesOnCell(:)
      integer, allocatable :: cellsOnCell(:,:)
      integer, allocatable :: edgesOnCell(:,:)
      real(kind=RKIND), allocatable :: h_weight(:,:)
      real(kind=RKIND), allocatable :: v_weight_upper(:,:)
      real(kind=RKIND), allocatable :: v_weight_lower(:,:)
      real(kind=RKIND), allocatable :: diag(:,:)
      real(kind=RKIND), allocatable :: volume(:,:)
      logical, allocatable :: active(:,:)
   contains
      procedure :: matvec => elliptic_operator_matvec
      procedure :: destroy => elliptic_operator_destroy
   end type mpas_elliptic_operator_type

contains

!***********************************************************************
!
!  routine elliptic_operator_init
!
!> \brief Initialize an electrostatic elliptic operator
!
!-----------------------------------------------------------------------
   subroutine elliptic_operator_init(op, nCells, nEdges, nVertLevels, maxEdges, &
                                    areaCell, dcEdge, dvEdge, nEdgesOnCell, &
                                    cellsOnCell, edgesOnCell, dz, epsilon0)

      type(mpas_elliptic_operator_type), intent(out) :: op
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      integer, intent(in) :: maxEdges
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: edgesOnCell
      real(kind=RKIND), dimension(nCells), intent(in) :: areaCell
      real(kind=RKIND), dimension(nEdges), intent(in) :: dcEdge
      real(kind=RKIND), dimension(nEdges), intent(in) :: dvEdge
      real(kind=RKIND), dimension(nVertLevels), intent(in) :: dz
      real(kind=RKIND), intent(in) :: epsilon0

      integer :: iCell, iEdge, k
      real(kind=RKIND) :: dzFace
      real(kind=RKIND), allocatable :: h_flat(:,:), v_up_flat(:,:), v_lo_flat(:,:)
      real(kind=RKIND), allocatable :: diag_extra_flat(:,:), volume_flat(:,:)
      logical, allocatable :: active_flat(:,:)

      allocate(h_flat(nVertLevels,nEdges))
      allocate(v_up_flat(nVertLevels,nCells))
      allocate(v_lo_flat(nVertLevels,nCells))
      allocate(diag_extra_flat(nVertLevels,nCells))
      allocate(volume_flat(nVertLevels,nCells))
      allocate(active_flat(nVertLevels,nCells))

      h_flat = 0.0_RKIND
      v_up_flat = 0.0_RKIND
      v_lo_flat = 0.0_RKIND
      diag_extra_flat = 0.0_RKIND
      volume_flat = 0.0_RKIND
      active_flat = .true.

      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            h_flat(k,iEdge) = epsilon0 * dvEdge(iEdge) * dz(k) / dcEdge(iEdge)
         end do
      end do

      do iCell = 1, nCells
         do k = 1, nVertLevels
            volume_flat(k,iCell) = areaCell(iCell) * dz(k)

            if (k < nVertLevels) then
               dzFace = 0.5_RKIND * (dz(k) + dz(k+1))
               v_up_flat(k,iCell) = epsilon0 * areaCell(iCell) / dzFace
            end if

            if (k > 1) then
               dzFace = 0.5_RKIND * (dz(k-1) + dz(k))
               v_lo_flat(k,iCell) = epsilon0 * areaCell(iCell) / dzFace
            else
               v_lo_flat(k,iCell) = epsilon0 * areaCell(iCell) / (0.5_RKIND * dz(k))
            end if
         end do
      end do

      call elliptic_operator_init_from_weights(op, nCells, nEdges, nVertLevels, maxEdges, &
                                               nEdgesOnCell, cellsOnCell, edgesOnCell, &
                                               h_flat, v_up_flat, v_lo_flat, diag_extra_flat, &
                                               volume_flat, active_flat)

      deallocate(h_flat, v_up_flat, v_lo_flat, diag_extra_flat, volume_flat, active_flat)

   end subroutine elliptic_operator_init

!***********************************************************************
!
!  routine elliptic_operator_init_from_weights
!
!> \brief Initialize an electrostatic elliptic operator from precomputed weights
!
!-----------------------------------------------------------------------
   subroutine elliptic_operator_init_from_weights(op, nCells, nEdges, nVertLevels, maxEdges, &
                                                  nEdgesOnCell, cellsOnCell, edgesOnCell, &
                                                  h_weight, v_weight_upper, v_weight_lower, &
                                                  diag_extra, volume, active)

      type(mpas_elliptic_operator_type), intent(out) :: op
      integer, intent(in) :: nCells
      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      integer, intent(in) :: maxEdges
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
      integer, dimension(maxEdges,nCells), intent(in) :: edgesOnCell
      real(kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: h_weight
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_upper
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: v_weight_lower
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: diag_extra
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: volume
      logical, dimension(nVertLevels,nCells), intent(in) :: active

      integer :: iCell, iEdge, iEdgeOnCell, jCell, k

      op%nCells = nCells
      op%nEdges = nEdges
      op%nVertLevels = nVertLevels
      op%maxEdges = maxEdges

      allocate(op%nEdgesOnCell(nCells))
      allocate(op%cellsOnCell(maxEdges,nCells))
      allocate(op%edgesOnCell(maxEdges,nCells))
      allocate(op%h_weight(nVertLevels,nEdges))
      allocate(op%v_weight_upper(nVertLevels,nCells))
      allocate(op%v_weight_lower(nVertLevels,nCells))
      allocate(op%diag(nVertLevels,nCells))
      allocate(op%volume(nVertLevels,nCells))
      allocate(op%active(nVertLevels,nCells))

      op%nEdgesOnCell = nEdgesOnCell
      op%cellsOnCell = cellsOnCell
      op%edgesOnCell = edgesOnCell
      op%h_weight = h_weight
      op%v_weight_upper = v_weight_upper
      op%v_weight_lower = v_weight_lower
      op%volume = volume
      op%active = active

      do iCell = 1, nCells
         do k = 1, nVertLevels
            if (.not. op%active(k,iCell)) then
               op%diag(k,iCell) = 1.0_RKIND
               cycle
            end if

            op%diag(k,iCell) = op%v_weight_upper(k,iCell) + op%v_weight_lower(k,iCell) + &
                               diag_extra(k,iCell)

            do iEdgeOnCell = 1, op%nEdgesOnCell(iCell)
               iEdge = op%edgesOnCell(iEdgeOnCell,iCell)
               jCell = op%cellsOnCell(iEdgeOnCell,iCell)
               if (iEdge > 0 .and. iEdge <= nEdges .and. &
                   jCell > 0 .and. jCell <= nCells) then
                  op%diag(k,iCell) = op%diag(k,iCell) + op%h_weight(k,iEdge)
               end if
            end do
         end do
      end do

   end subroutine elliptic_operator_init_from_weights

!***********************************************************************
!
!  routine elliptic_operator_matvec
!
!> \brief Apply the electrostatic elliptic operator
!
!-----------------------------------------------------------------------
   subroutine elliptic_operator_matvec(op, x, y)

      class(mpas_elliptic_operator_type), intent(in) :: op
      real(kind=RKIND), dimension(op%nVertLevels,op%nCells), intent(in) :: x
      real(kind=RKIND), dimension(op%nVertLevels,op%nCells), intent(out) :: y

      integer :: iCell, iEdge, iEdgeOnCell, jCell, k

      y = 0.0_RKIND

      do iCell = 1, op%nCells
         do k = 1, op%nVertLevels
            y(k,iCell) = op%diag(k,iCell) * x(k,iCell)
         end do

         do iEdgeOnCell = 1, op%nEdgesOnCell(iCell)
            iEdge = op%edgesOnCell(iEdgeOnCell,iCell)
            jCell = op%cellsOnCell(iEdgeOnCell,iCell)
            if (iEdge > 0 .and. iEdge <= op%nEdges .and. &
                jCell > 0 .and. jCell <= op%nCells) then
               do k = 1, op%nVertLevels
                  if (op%active(k,iCell) .and. op%active(k,jCell)) &
                     y(k,iCell) = y(k,iCell) - op%h_weight(k,iEdge) * x(k,jCell)
               end do
            end if
         end do

         do k = 1, op%nVertLevels - 1
            if (op%active(k,iCell) .and. op%active(k+1,iCell)) &
               y(k,iCell) = y(k,iCell) - op%v_weight_upper(k,iCell) * x(k+1,iCell)
         end do

         do k = 2, op%nVertLevels
            if (op%active(k,iCell) .and. op%active(k-1,iCell)) &
               y(k,iCell) = y(k,iCell) - op%v_weight_lower(k,iCell) * x(k-1,iCell)
         end do
      end do

   end subroutine elliptic_operator_matvec

!***********************************************************************
!
!  routine elliptic_operator_destroy
!
!> \brief Release storage owned by an electrostatic elliptic operator
!
!-----------------------------------------------------------------------
   subroutine elliptic_operator_destroy(op)

      class(mpas_elliptic_operator_type), intent(inout) :: op

      if (allocated(op%nEdgesOnCell)) deallocate(op%nEdgesOnCell)
      if (allocated(op%cellsOnCell)) deallocate(op%cellsOnCell)
      if (allocated(op%edgesOnCell)) deallocate(op%edgesOnCell)
      if (allocated(op%h_weight)) deallocate(op%h_weight)
      if (allocated(op%v_weight_upper)) deallocate(op%v_weight_upper)
      if (allocated(op%v_weight_lower)) deallocate(op%v_weight_lower)
      if (allocated(op%diag)) deallocate(op%diag)
      if (allocated(op%volume)) deallocate(op%volume)
      if (allocated(op%active)) deallocate(op%active)

      op%nCells = 0
      op%nEdges = 0
      op%nVertLevels = 0
      op%maxEdges = 0

   end subroutine elliptic_operator_destroy

end module mpas_elliptic_operator
