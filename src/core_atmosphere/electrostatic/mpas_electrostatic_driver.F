! Copyright (c) 2026,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_electrostatic_driver
!
!> \brief Diagnostic electrostatic Poisson driver for MPAS-Atmosphere
!
!-----------------------------------------------------------------------
module mpas_electrostatic_driver

   use mpas_derived_types
   use mpas_kind_types, only : RKIND, StrKIND
   use mpas_pool_routines
   use mpas_log, only : mpas_log_write
   use mpas_atm_halos, only : exchange_halo_group
   use mpas_vector_reconstruction, only : mpas_reconstruct
   use mpas_elliptic_operator, only : mpas_elliptic_operator_type, elliptic_operator_init
   use mpas_pcg_solver, only : mpas_pcg_solver_type, pcg_init, PC_NONE, PC_JACOBI
   use mpas_electrostatic_bcs, only : electrostatic_apply_ground_rhs, &
                                      electrostatic_compute_vertical_E
   use mpas_electrostatic_source, only : electrostatic_fill_mms_cart_source

   implicit none
   private
   save

   public :: mpas_electrostatic_init
   public :: mpas_electrostatic_solve
   public :: mpas_electrostatic_finalize

   real(kind=RKIND), parameter :: epsilon0 = 8.8541878128e-12_RKIND

   type electrostatic_block_state
      type(mpas_elliptic_operator_type) :: op
      type(mpas_pcg_solver_type) :: solver
      real(kind=RKIND), allocatable :: rhs(:,:)
      real(kind=RKIND), allocatable :: zMid(:,:)
      real(kind=RKIND), allocatable :: E_vertical(:,:)
      real(kind=RKIND), allocatable :: E_x(:,:)
      real(kind=RKIND), allocatable :: E_y(:,:)
      real(kind=RKIND), allocatable :: E_z(:,:)
      real(kind=RKIND), allocatable :: E_zonal(:,:)
      real(kind=RKIND), allocatable :: E_meridional(:,:)
      real(kind=RKIND) :: x_period = 1.0_RKIND
      real(kind=RKIND) :: y_period = 1.0_RKIND
      real(kind=RKIND) :: z_top = 1.0_RKIND
      type(electrostatic_block_state), pointer :: next => null()
   end type electrostatic_block_state

   type(electrostatic_block_state), pointer :: block_states => null()
   logical :: electrostatic_initialized = .false.

contains

!***********************************************************************
!
!  routine mpas_electrostatic_init
!
!> \brief Initialize diagnostic electrostatic workspaces for each block
!
!-----------------------------------------------------------------------
   subroutine mpas_electrostatic_init(domain)

      type(domain_type), intent(inout) :: domain

      type(block_type), pointer :: block_ptr
      type(electrostatic_block_state), pointer :: new_state, tail

      if (electrostatic_initialized) then
         call mpas_electrostatic_finalize(domain)
      end if

      nullify(block_states)
      nullify(tail)

      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         allocate(new_state)
         nullify(new_state % next)

         call init_block_state(block_ptr, new_state)

         if (.not. associated(block_states)) then
            block_states => new_state
         else
            tail % next => new_state
         end if
         tail => new_state

         block_ptr => block_ptr % next
      end do

      electrostatic_initialized = .true.

   end subroutine mpas_electrostatic_init

!***********************************************************************
!
!  routine mpas_electrostatic_solve
!
!> \brief Solve diagnostic electrostatic Poisson equation and reconstruct E
!
!-----------------------------------------------------------------------
   subroutine mpas_electrostatic_solve(domain)

      type(domain_type), intent(inout) :: domain

      type(block_type), pointer :: block_ptr
      type(electrostatic_block_state), pointer :: state_ptr

      if (.not. electrostatic_initialized) then
         call mpas_electrostatic_init(domain)
      end if

      block_ptr => domain % blocklist
      state_ptr => block_states
      do while (associated(block_ptr) .and. associated(state_ptr))
         call solve_block_potential(block_ptr, state_ptr)
         block_ptr => block_ptr % next
         state_ptr => state_ptr % next
      end do

      call exchange_halo_group(domain, 'electrostatic:phi')

      block_ptr => domain % blocklist
      state_ptr => block_states
      do while (associated(block_ptr) .and. associated(state_ptr))
         call postprocess_block_E(block_ptr, state_ptr)
         block_ptr => block_ptr % next
         state_ptr => state_ptr % next
      end do

   end subroutine mpas_electrostatic_solve

!***********************************************************************
!
!  routine mpas_electrostatic_finalize
!
!> \brief Release diagnostic electrostatic workspaces
!
!-----------------------------------------------------------------------
   subroutine mpas_electrostatic_finalize(domain)

      type(domain_type), intent(inout) :: domain

      type(electrostatic_block_state), pointer :: state_ptr, next_state

      state_ptr => block_states
      do while (associated(state_ptr))
         next_state => state_ptr % next

         call state_ptr % op % destroy()
         call state_ptr % solver % destroy()
         if (allocated(state_ptr % rhs)) deallocate(state_ptr % rhs)
         if (allocated(state_ptr % zMid)) deallocate(state_ptr % zMid)
         if (allocated(state_ptr % E_vertical)) deallocate(state_ptr % E_vertical)
         if (allocated(state_ptr % E_x)) deallocate(state_ptr % E_x)
         if (allocated(state_ptr % E_y)) deallocate(state_ptr % E_y)
         if (allocated(state_ptr % E_z)) deallocate(state_ptr % E_z)
         if (allocated(state_ptr % E_zonal)) deallocate(state_ptr % E_zonal)
         if (allocated(state_ptr % E_meridional)) deallocate(state_ptr % E_meridional)

         deallocate(state_ptr)
         state_ptr => next_state
      end do

      nullify(block_states)
      electrostatic_initialized = .false.

   end subroutine mpas_electrostatic_finalize

!***********************************************************************
!
!  routine init_block_state
!
!> \brief Initialize one block's electrostatic operator and PCG workspace
!
!-----------------------------------------------------------------------
   subroutine init_block_state(block_ptr, state_ptr)

      type(block_type), pointer, intent(inout) :: block_ptr
      type(electrostatic_block_state), pointer, intent(inout) :: state_ptr

      type(mpas_pool_type), pointer :: mesh_pool, es_pool
      integer, pointer :: nCells_ptr, nEdges_ptr, nVertLevels_ptr, maxEdges_ptr
      integer, pointer :: max_iter
      integer :: nCells, nEdges, nVertLevels, maxEdges, preconditioner
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell
      real(kind=RKIND), dimension(:), pointer :: areaCell, dcEdge, dvEdge
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell
      real(kind=RKIND), dimension(:,:), pointer :: zgrid
      real(kind=RKIND), dimension(:,:), pointer :: phi, rho_charge, E_normal
      real(kind=RKIND), dimension(:,:,:), pointer :: E_vector
      real(kind=RKIND), pointer :: tol
      real(kind=RKIND), allocatable :: dz(:)
      character(len=StrKIND), pointer :: pc_name

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
      call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)

      call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nEdges', nEdges_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'maxEdges', maxEdges_ptr)
      nCells = nCells_ptr
      nEdges = nEdges_ptr
      nVertLevels = nVertLevels_ptr
      maxEdges = maxEdges_ptr

      call mpas_pool_get_array(mesh_pool, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh_pool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh_pool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh_pool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh_pool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh_pool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh_pool, 'xCell', xCell)
      call mpas_pool_get_array(mesh_pool, 'yCell', yCell)
      call mpas_pool_get_array(mesh_pool, 'zgrid', zgrid)

      call mpas_pool_get_config(block_ptr % configs, 'config_poisson_preconditioner', pc_name)
      call mpas_pool_get_config(block_ptr % configs, 'config_poisson_tol', tol)
      call mpas_pool_get_config(block_ptr % configs, 'config_poisson_max_iter', max_iter)

      allocate(state_ptr % rhs(nVertLevels,nCells))
      allocate(state_ptr % zMid(nVertLevels,nCells))
      allocate(state_ptr % E_vertical(nVertLevels,nCells))
      allocate(state_ptr % E_x(nVertLevels,nCells))
      allocate(state_ptr % E_y(nVertLevels,nCells))
      allocate(state_ptr % E_z(nVertLevels,nCells))
      allocate(state_ptr % E_zonal(nVertLevels,nCells))
      allocate(state_ptr % E_meridional(nVertLevels,nCells))
      allocate(dz(nVertLevels))

      call fill_z_geometry(nCells, nVertLevels, zgrid, state_ptr % zMid, dz)
      call infer_mms_periods(nCells, nVertLevels, xCell, yCell, zgrid, &
                             state_ptr % x_period, state_ptr % y_period, state_ptr % z_top)

      preconditioner = preconditioner_from_name(pc_name)
      call elliptic_operator_init(state_ptr % op, nCells, nEdges, nVertLevels, maxEdges, &
                                  areaCell, dcEdge, dvEdge, nEdgesOnCell, cellsOnCell, &
                                  edgesOnCell, dz, epsilon0)
      call pcg_init(state_ptr % solver, nCells, nVertLevels, preconditioner, max(1, max_iter), tol)

      state_ptr % rhs = 0.0_RKIND
      state_ptr % E_vertical = 0.0_RKIND
      state_ptr % E_x = 0.0_RKIND
      state_ptr % E_y = 0.0_RKIND
      state_ptr % E_z = 0.0_RKIND
      state_ptr % E_zonal = 0.0_RKIND
      state_ptr % E_meridional = 0.0_RKIND

      call mpas_pool_get_array(es_pool, 'phi', phi)
      call mpas_pool_get_array(es_pool, 'rho_charge', rho_charge)
      call mpas_pool_get_array(es_pool, 'E_normal', E_normal)
      call mpas_pool_get_array(es_pool, 'E_vector', E_vector)
      phi = 0.0_RKIND
      rho_charge = 0.0_RKIND
      E_normal = 0.0_RKIND
      E_vector = 0.0_RKIND
      call set_solver_diagnostics(es_pool, 0, 0.0_RKIND, 0.0_RKIND)

      deallocate(dz)

   end subroutine init_block_state

!***********************************************************************
!
!  routine solve_block_potential
!
!> \brief Fill source terms and solve A phi = rhs for one block
!
!-----------------------------------------------------------------------
   subroutine solve_block_potential(block_ptr, state_ptr)

      type(block_type), pointer, intent(inout) :: block_ptr
      type(electrostatic_block_state), pointer, intent(inout) :: state_ptr

      type(mpas_pool_type), pointer :: mesh_pool, es_pool
      integer, pointer :: nCells_ptr, nVertLevels_ptr
      integer :: nCells, nVertLevels, iter_count, status
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell
      real(kind=RKIND), dimension(:,:), pointer :: zgrid
      real(kind=RKIND), dimension(:,:), pointer :: phi, rho_charge
      real(kind=RKIND), pointer :: phi_ground
      real(kind=RKIND) :: r0, rfinal
      character(len=StrKIND), pointer :: source_name

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
      call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)

      call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels_ptr)
      nCells = nCells_ptr
      nVertLevels = nVertLevels_ptr

      call mpas_pool_get_array(mesh_pool, 'xCell', xCell)
      call mpas_pool_get_array(mesh_pool, 'yCell', yCell)
      call mpas_pool_get_array(mesh_pool, 'zgrid', zgrid)
      call mpas_pool_get_array(es_pool, 'phi', phi)
      call mpas_pool_get_array(es_pool, 'rho_charge', rho_charge)

      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_source', source_name)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_bc_ground', phi_ground)

      call fill_z_mid(nCells, nVertLevels, zgrid, state_ptr % zMid)

      if (trim(source_name) /= 'external') then
         call electrostatic_fill_mms_cart_source(nCells, nVertLevels, xCell, yCell, state_ptr % zMid, &
                                                 epsilon0, state_ptr % x_period, state_ptr % y_period, &
                                                 state_ptr % z_top, rho_charge)
      end if

      state_ptr % rhs = rho_charge * state_ptr % op % volume

      if (abs(phi_ground) > 0.0_RKIND) then
         call electrostatic_apply_ground_rhs(nCells, nVertLevels, phi_ground, &
                                             state_ptr % op % v_weight_lower, state_ptr % rhs)
      end if

      call state_ptr % solver % solve(state_ptr % op, state_ptr % rhs, phi, &
                                      iter_count, r0, rfinal, status)

      call set_solver_diagnostics(es_pool, iter_count, r0, rfinal)

      if (status /= 0) then
         call mpas_log_write('Electrostatic PCG returned status $i after $i iterations; final residual = $r', &
                             messageType=MPAS_LOG_WARN, intArgs=(/status, iter_count/), realArgs=(/rfinal/))
      end if

   end subroutine solve_block_potential

!***********************************************************************
!
!  routine postprocess_block_E
!
!> \brief Compute edge-normal and cell-centered electric fields
!
!-----------------------------------------------------------------------
   subroutine postprocess_block_E(block_ptr, state_ptr)

      type(block_type), pointer, intent(inout) :: block_ptr
      type(electrostatic_block_state), pointer, intent(inout) :: state_ptr

      type(mpas_pool_type), pointer :: mesh_pool, es_pool
      integer, pointer :: nCells_ptr, nEdges_ptr, nVertLevels_ptr
      integer :: iCell, k, nCells, nEdges, nVertLevels
      integer, dimension(:,:), pointer :: cellsOnEdge
      real(kind=RKIND), dimension(:), pointer :: dcEdge
      real(kind=RKIND), dimension(:,:), pointer :: localVerticalUnitVectors
      real(kind=RKIND), dimension(:,:), pointer :: zgrid
      real(kind=RKIND), dimension(:,:), pointer :: phi, E_normal
      real(kind=RKIND), dimension(:,:,:), pointer :: E_vector
      real(kind=RKIND), pointer :: phi_ground

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
      call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)

      call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nEdges', nEdges_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels_ptr)
      nCells = nCells_ptr
      nEdges = nEdges_ptr
      nVertLevels = nVertLevels_ptr

      call mpas_pool_get_array(mesh_pool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh_pool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh_pool, 'localVerticalUnitVectors', localVerticalUnitVectors)
      call mpas_pool_get_array(mesh_pool, 'zgrid', zgrid)
      call mpas_pool_get_array(es_pool, 'phi', phi)
      call mpas_pool_get_array(es_pool, 'E_normal', E_normal)
      call mpas_pool_get_array(es_pool, 'E_vector', E_vector)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_bc_ground', phi_ground)

      call compute_edge_normal_E(nEdges, nVertLevels, cellsOnEdge, dcEdge, phi, E_normal)

      call mpas_reconstruct(mesh_pool, E_normal, state_ptr % E_x, state_ptr % E_y, state_ptr % E_z, &
                            state_ptr % E_zonal, state_ptr % E_meridional, includeHalos=.true.)

      call fill_z_mid(nCells, nVertLevels, zgrid, state_ptr % zMid)
      call electrostatic_compute_vertical_E(nCells, nVertLevels, phi, state_ptr % zMid, &
                                            phi_ground, state_ptr % E_vertical)

      do iCell = 1, nCells
         do k = 1, nVertLevels
            E_vector(1,k,iCell) = state_ptr % E_x(k,iCell) + &
                                  state_ptr % E_vertical(k,iCell) * localVerticalUnitVectors(1,iCell)
            E_vector(2,k,iCell) = state_ptr % E_y(k,iCell) + &
                                  state_ptr % E_vertical(k,iCell) * localVerticalUnitVectors(2,iCell)
            E_vector(3,k,iCell) = state_ptr % E_z(k,iCell) + &
                                  state_ptr % E_vertical(k,iCell) * localVerticalUnitVectors(3,iCell)
         end do
      end do

   end subroutine postprocess_block_E

!***********************************************************************
!
!  routine compute_edge_normal_E
!
!> \brief Compute lateral edge-normal E from cell-centered potential
!
!-----------------------------------------------------------------------
   subroutine compute_edge_normal_E(nEdges, nVertLevels, cellsOnEdge, dcEdge, phi, E_normal)

      integer, intent(in) :: nEdges
      integer, intent(in) :: nVertLevels
      integer, dimension(2,nEdges), intent(in) :: cellsOnEdge
      real(kind=RKIND), dimension(nEdges), intent(in) :: dcEdge
      real(kind=RKIND), dimension(nVertLevels,*), intent(in) :: phi
      real(kind=RKIND), dimension(nVertLevels,nEdges), intent(out) :: E_normal

      integer :: iEdge, cell1, cell2, k

      E_normal = 0.0_RKIND

      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 > 0 .and. cell2 > 0 .and. dcEdge(iEdge) > 0.0_RKIND) then
            do k = 1, nVertLevels
               E_normal(k,iEdge) = -(phi(k,cell2) - phi(k,cell1)) / dcEdge(iEdge)
            end do
         end if
      end do

   end subroutine compute_edge_normal_E

!***********************************************************************
!
!  routine fill_z_geometry
!
!> \brief Compute cell-centered z and representative layer thicknesses
!
!-----------------------------------------------------------------------
   subroutine fill_z_geometry(nCells, nVertLevels, zgrid, zMid, dz)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels+1,nCells), intent(in) :: zgrid
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: zMid
      real(kind=RKIND), dimension(nVertLevels), intent(out) :: dz

      integer :: k

      call fill_z_mid(nCells, nVertLevels, zgrid, zMid)

      do k = 1, nVertLevels
         dz(k) = sum(zgrid(k+1,1:nCells) - zgrid(k,1:nCells)) / real(nCells, kind=RKIND)
         if (dz(k) <= 0.0_RKIND) then
            dz(k) = 1.0_RKIND
         end if
      end do

   end subroutine fill_z_geometry

!***********************************************************************
!
!  routine fill_z_mid
!
!> \brief Compute cell-centered z from interface heights
!
!-----------------------------------------------------------------------
   subroutine fill_z_mid(nCells, nVertLevels, zgrid, zMid)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels+1,nCells), intent(in) :: zgrid
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(out) :: zMid

      integer :: iCell, k

      do iCell = 1, nCells
         do k = 1, nVertLevels
            zMid(k,iCell) = 0.5_RKIND * (zgrid(k,iCell) + zgrid(k+1,iCell))
         end do
      end do

   end subroutine fill_z_mid

!***********************************************************************
!
!  routine infer_mms_periods
!
!> \brief Infer Cartesian MMS length scales from available mesh geometry
!
!-----------------------------------------------------------------------
   subroutine infer_mms_periods(nCells, nVertLevels, xCell, yCell, zgrid, x_period, y_period, z_top)

      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nCells), intent(in) :: xCell
      real(kind=RKIND), dimension(nCells), intent(in) :: yCell
      real(kind=RKIND), dimension(nVertLevels+1,nCells), intent(in) :: zgrid
      real(kind=RKIND), intent(out) :: x_period
      real(kind=RKIND), intent(out) :: y_period
      real(kind=RKIND), intent(out) :: z_top

      x_period = maxval(xCell(1:nCells)) - minval(xCell(1:nCells))
      y_period = maxval(yCell(1:nCells)) - minval(yCell(1:nCells))
      z_top = maxval(zgrid(1:nVertLevels+1,1:nCells))

      if (x_period <= 0.0_RKIND) x_period = 1.0_RKIND
      if (y_period <= 0.0_RKIND) y_period = 1.0_RKIND
      if (z_top <= 0.0_RKIND) z_top = 1.0_RKIND

   end subroutine infer_mms_periods

!***********************************************************************
!
!  routine preconditioner_from_name
!
!> \brief Map namelist preconditioner strings to PCG constants
!
!-----------------------------------------------------------------------
   integer function preconditioner_from_name(pc_name)

      character(len=*), intent(in) :: pc_name

      select case (trim(pc_name))
      case ('none')
         preconditioner_from_name = PC_NONE
      case default
         preconditioner_from_name = PC_JACOBI
      end select

   end function preconditioner_from_name

!***********************************************************************
!
!  routine set_solver_diagnostics
!
!> \brief Store scalar solver diagnostics in electrostatic fields
!
!-----------------------------------------------------------------------
   subroutine set_solver_diagnostics(es_pool, iter_count, r0, rfinal)

      type(mpas_pool_type), pointer, intent(in) :: es_pool
      integer, intent(in) :: iter_count
      real(kind=RKIND), intent(in) :: r0
      real(kind=RKIND), intent(in) :: rfinal

      type(field0DInteger), pointer :: iter_field
      type(field0DReal), pointer :: r0_field, rf_field

      call mpas_pool_get_field(es_pool, 'cg_iter_count', iter_field)
      call mpas_pool_get_field(es_pool, 'cg_residual_initial', r0_field)
      call mpas_pool_get_field(es_pool, 'cg_residual_final', rf_field)

      iter_field % scalar = iter_count
      r0_field % scalar = r0
      rf_field % scalar = rfinal

   end subroutine set_solver_diagnostics

end module mpas_electrostatic_driver
