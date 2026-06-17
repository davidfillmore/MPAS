! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_pcg_solver
!
!> \brief Preconditioned conjugate-gradient solver kernels (serial reference)
!> \details
!>  Single-partition reference implementation: the dot products are local sums
!>  with no MPI reduction and the matvec performs no halo exchange, so it is
!>  valid only when the entire problem lives on one partition
!>  (nCells == nCellsSolve), i.e. the unit tests. The MPI-parallel solve used
!>  in production is solve_domain_potential in mpas_electrostatic_driver, which
!>  halo-exchanges the iterate before each matvec and globally reduces every
!>  inner product over owned cells only.
!
!-----------------------------------------------------------------------
module mpas_pcg_solver

   use mpas_kind_types, only : RKIND
   use mpas_elliptic_operator, only : mpas_elliptic_operator_type

   implicit none
   private
   save

   integer, parameter, public :: PC_NONE = 0
   integer, parameter, public :: PC_JACOBI = 1

   public :: mpas_pcg_solver_type, pcg_init

   type :: mpas_pcg_solver_type
      integer :: nCells = 0
      integer :: nVertLevels = 0
      integer :: preconditioner = PC_NONE
      integer :: max_iter = 0
      real(kind=RKIND) :: tol = 0.0_RKIND
      real(kind=RKIND), allocatable :: r(:,:)
      real(kind=RKIND), allocatable :: z(:,:)
      real(kind=RKIND), allocatable :: p(:,:)
      real(kind=RKIND), allocatable :: Ap(:,:)
   contains
      procedure :: solve => pcg_solve
      procedure :: destroy => pcg_destroy
   end type mpas_pcg_solver_type

contains

!***********************************************************************
!
!  routine pcg_init
!
!> \brief Initialize a PCG solver workspace
!
!-----------------------------------------------------------------------
   subroutine pcg_init(solver, nCells, nVertLevels, preconditioner, max_iter, tol)

      type(mpas_pcg_solver_type), intent(out) :: solver
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      integer, intent(in) :: preconditioner
      integer, intent(in) :: max_iter
      real(kind=RKIND), intent(in) :: tol

      solver%nCells = nCells
      solver%nVertLevels = nVertLevels
      solver%preconditioner = preconditioner
      solver%max_iter = max_iter
      solver%tol = tol

      allocate(solver%r(nVertLevels,nCells))
      allocate(solver%z(nVertLevels,nCells))
      allocate(solver%p(nVertLevels,nCells))
      allocate(solver%Ap(nVertLevels,nCells))

      solver%r = 0.0_RKIND
      solver%z = 0.0_RKIND
      solver%p = 0.0_RKIND
      solver%Ap = 0.0_RKIND

   end subroutine pcg_init

!***********************************************************************
!
!  routine pcg_solve
!
!> \brief Solve A phi = rhs with preconditioned conjugate gradients (serial)
!> \details
!>  Single-partition reference path (no halo exchange, no MPI reduction); see
!>  the module header. status = 0 on convergence, 1 on max iteration, 2 on
!>  non-positive curvature, and 3 on invalid preconditioner data.
!
!-----------------------------------------------------------------------
   subroutine pcg_solve(solver, op, rhs, phi, iter_count, r0, rfinal, status)

      class(mpas_pcg_solver_type), intent(inout) :: solver
      type(mpas_elliptic_operator_type), intent(in) :: op
      real(kind=RKIND), dimension(solver%nVertLevels,solver%nCells), intent(in) :: rhs
      real(kind=RKIND), dimension(solver%nVertLevels,solver%nCells), intent(inout) :: phi
      integer, intent(out) :: iter_count
      real(kind=RKIND), intent(out) :: r0
      real(kind=RKIND), intent(out) :: rfinal
      integer, intent(out) :: status

      integer :: iter
      real(kind=RKIND) :: alpha, beta, denom, rhs_norm, rz_new, rz_old

      iter_count = 0
      status = 1

      call op%matvec(phi, solver%Ap)
      solver%r = rhs - solver%Ap

      rhs_norm = sqrt(dot2(rhs, rhs))
      r0 = residual_norm(solver%r, rhs_norm)
      rfinal = r0

      if (rfinal <= solver%tol) then
         status = 0
         return
      end if

      if (.not. apply_preconditioner(solver, op, status)) then
         return
      end if

      rz_old = dot2(solver%r, solver%z)
      if (rz_old <= 0.0_RKIND) then
         status = 2
         return
      end if

      solver%p = solver%z

      do iter = 1, solver%max_iter
         call op%matvec(solver%p, solver%Ap)
         denom = dot2(solver%p, solver%Ap)
         if (denom <= 0.0_RKIND) then
            iter_count = iter - 1
            status = 2
            return
         end if

         alpha = rz_old / denom
         phi = phi + alpha * solver%p
         solver%r = solver%r - alpha * solver%Ap

         rfinal = residual_norm(solver%r, rhs_norm)
         iter_count = iter
         if (rfinal <= solver%tol) then
            status = 0
            return
         end if

         if (.not. apply_preconditioner(solver, op, status)) then
            return
         end if

         rz_new = dot2(solver%r, solver%z)
         if (rz_new <= 0.0_RKIND) then
            status = 2
            return
         end if

         beta = rz_new / rz_old
         solver%p = solver%z + beta * solver%p
         rz_old = rz_new
      end do

      status = 1

   end subroutine pcg_solve

!***********************************************************************
!
!  routine apply_preconditioner
!
!> \brief Apply the configured PCG preconditioner
!
!-----------------------------------------------------------------------
   logical function apply_preconditioner(solver, op, status)

      class(mpas_pcg_solver_type), intent(inout) :: solver
      type(mpas_elliptic_operator_type), intent(in) :: op
      integer, intent(out) :: status

      apply_preconditioner = .false.

      select case (solver%preconditioner)
      case (PC_NONE)
         solver%z = solver%r
      case (PC_JACOBI)
         if (any(op%diag <= 0.0_RKIND)) then
            status = 3
            return
         end if
         solver%z = solver%r / op%diag
      case default
         status = 3
         return
      end select

      status = 1
      apply_preconditioner = .true.

   end function apply_preconditioner

!***********************************************************************
!
!  routine residual_norm
!
!> \brief Return a relative residual norm, falling back to absolute norm
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function residual_norm(r, rhs_norm)

      real(kind=RKIND), dimension(:,:), intent(in) :: r
      real(kind=RKIND), intent(in) :: rhs_norm
      real(kind=RKIND) :: r_norm

      r_norm = sqrt(dot2(r, r))
      if (rhs_norm > 0.0_RKIND) then
         residual_norm = r_norm / rhs_norm
      else
         residual_norm = r_norm
      end if

   end function residual_norm

!***********************************************************************
!
!  routine dot2
!
!> \brief Dot product for rank-2 solver work arrays
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function dot2(a, b)

      real(kind=RKIND), dimension(:,:), intent(in) :: a
      real(kind=RKIND), dimension(:,:), intent(in) :: b

      dot2 = sum(a * b)

   end function dot2

!***********************************************************************
!
!  routine pcg_destroy
!
!> \brief Release PCG solver workspace
!
!-----------------------------------------------------------------------
   subroutine pcg_destroy(solver)

      class(mpas_pcg_solver_type), intent(inout) :: solver

      if (allocated(solver%r)) deallocate(solver%r)
      if (allocated(solver%z)) deallocate(solver%z)
      if (allocated(solver%p)) deallocate(solver%p)
      if (allocated(solver%Ap)) deallocate(solver%Ap)

      solver%nCells = 0
      solver%nVertLevels = 0
      solver%preconditioner = PC_NONE
      solver%max_iter = 0
      solver%tol = 0.0_RKIND

   end subroutine pcg_destroy

end module mpas_pcg_solver
