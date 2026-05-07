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

      integer :: iCell, iEdge, iEdgeOnCell, jCell, k
      real(kind=RKIND) :: dzFace

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

      op%nEdgesOnCell = nEdgesOnCell
      op%cellsOnCell = cellsOnCell
      op%edgesOnCell = edgesOnCell
      op%h_weight = 0.0_RKIND
      op%v_weight_upper = 0.0_RKIND
      op%v_weight_lower = 0.0_RKIND
      op%diag = 0.0_RKIND
      op%volume = 0.0_RKIND

      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            op%h_weight(k,iEdge) = epsilon0 * dvEdge(iEdge) * dz(k) / dcEdge(iEdge)
         end do
      end do

      do iCell = 1, nCells
         do k = 1, nVertLevels
            op%volume(k,iCell) = areaCell(iCell) * dz(k)

            if (k < nVertLevels) then
               dzFace = 0.5_RKIND * (dz(k) + dz(k+1))
               op%v_weight_upper(k,iCell) = epsilon0 * areaCell(iCell) / dzFace
            end if

            if (k > 1) then
               dzFace = 0.5_RKIND * (dz(k-1) + dz(k))
               op%v_weight_lower(k,iCell) = epsilon0 * areaCell(iCell) / dzFace
            else
               op%v_weight_lower(k,iCell) = epsilon0 * areaCell(iCell) / (0.5_RKIND * dz(k))
            end if
         end do
      end do

      do iCell = 1, nCells
         do k = 1, nVertLevels
            op%diag(k,iCell) = op%v_weight_upper(k,iCell) + op%v_weight_lower(k,iCell)

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

   end subroutine elliptic_operator_init

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

            do iEdgeOnCell = 1, op%nEdgesOnCell(iCell)
               iEdge = op%edgesOnCell(iEdgeOnCell,iCell)
               jCell = op%cellsOnCell(iEdgeOnCell,iCell)
               if (iEdge > 0 .and. iEdge <= op%nEdges .and. &
                   jCell > 0 .and. jCell <= op%nCells) then
                  y(k,iCell) = y(k,iCell) - op%h_weight(k,iEdge) * x(k,jCell)
               end if
            end do

            if (k < op%nVertLevels) then
               y(k,iCell) = y(k,iCell) - op%v_weight_upper(k,iCell) * x(k+1,iCell)
            end if

            if (k > 1) then
               y(k,iCell) = y(k,iCell) - op%v_weight_lower(k,iCell) * x(k-1,iCell)
            end if
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

      op%nCells = 0
      op%nEdges = 0
      op%nVertLevels = 0
      op%maxEdges = 0

   end subroutine elliptic_operator_destroy

end module mpas_elliptic_operator
