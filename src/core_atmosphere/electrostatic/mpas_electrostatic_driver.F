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
   use mpas_dmpar, only : mpas_dmpar_sum_real, mpas_dmpar_sum_real_array, &
                          mpas_dmpar_min_real, mpas_dmpar_max_real, &
                          mpas_dmpar_max_int
   use mpas_atm_halos, only : exchange_halo_group
   use mpas_vector_reconstruction, only : mpas_reconstruct
   use mpas_elliptic_operator, only : mpas_elliptic_operator_type, elliptic_operator_init_from_weights
   use mpas_pcg_solver, only : mpas_pcg_solver_type, pcg_init, PC_NONE, PC_JACOBI
   use mpas_electrostatic_bcs, only : electrostatic_apply_ground_rhs, &
                                      electrostatic_compute_edge_normal_E, &
                                      electrostatic_compute_vertical_E, &
                                      electrostatic_horizontal_metric_scale
   use mpas_electrostatic_source, only : electrostatic_fill_mms_cart_source, &
                                         electrostatic_fill_mms_cart_horizontal_source, &
                                         electrostatic_fill_mms_sphere_source, &
                                         electrostatic_fill_mms_sphere_horizontal_source, &
                                         electrostatic_fill_mms_terrain_source, &
                                         electrostatic_fill_dipole_source, &
                                         electrostatic_fill_tripole_source, &
                                         electrostatic_fill_gaussian_monopole_source, &
                                         electrostatic_fill_gaussian_dipole_y_source, &
                                         electrostatic_fill_gaussian_dipole_z_source, &
                                         electrostatic_fill_stub_source, &
                                         electrostatic_select_stub_indices
   use mpas_electrostatic_vertical_remap, only : electrostatic_build_poisson_grid, &
                                                 electrostatic_build_poisson_fv_weights, &
                                                 electrostatic_remap_rho_conservative, &
                                                 electrostatic_interp_poisson_to_mpas

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
      integer :: nCellsSolve = 0
      logical :: zgrid_mode = .false.
      real(kind=RKIND), allocatable :: rhs(:,:)
      real(kind=RKIND), allocatable :: x(:,:)
      real(kind=RKIND), allocatable :: zMid(:,:)
      real(kind=RKIND), allocatable :: zGround(:)
      real(kind=RKIND), allocatable :: zFaceP(:)
      real(kind=RKIND), allocatable :: zMidPC(:,:)
      real(kind=RKIND), allocatable :: airThickness(:,:)
      real(kind=RKIND), allocatable :: wallWeight(:,:)
      real(kind=RKIND), allocatable :: work(:,:)   ! per-solve scratch: rhs source buffer + back-remap target
      integer, allocatable :: kFirstActive(:)
      logical, allocatable :: active(:,:)
      real(kind=RKIND), allocatable :: E_vertical(:,:)
      real(kind=RKIND), allocatable :: E_x(:,:)
      real(kind=RKIND), allocatable :: E_y(:,:)
      real(kind=RKIND), allocatable :: E_z(:,:)
      real(kind=RKIND), allocatable :: E_zonal(:,:)
      real(kind=RKIND), allocatable :: E_meridional(:,:)
      real(kind=RKIND) :: x_period = 1.0_RKIND
      real(kind=RKIND) :: y_period = 1.0_RKIND
      real(kind=RKIND) :: z_top = 1.0_RKIND
      real(kind=RKIND) :: sphere_radius = 6371229.0_RKIND
      real(kind=RKIND) :: horizontal_length_scale = 1.0_RKIND
      real(kind=RKIND) :: source_x0 = 0.0_RKIND
      real(kind=RKIND) :: source_y0 = 0.0_RKIND
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
      character(len=StrKIND), pointer :: terrain_mode
      logical :: terrain_detected, zgrid_mode
      real(kind=RKIND) :: source_x0, source_y0
      real(kind=RKIND) :: mms_x, mms_y, z_top, z_bottom
      real(kind=RKIND), allocatable :: zFaceProfile(:)

      if (electrostatic_initialized) then
         call mpas_electrostatic_finalize(domain)
      end if

      nullify(block_states)
      nullify(tail)

      call compute_domain_synthetic_source_center(domain, source_x0, source_y0)
      call compute_domain_grid_extents(domain, mms_x, mms_y, z_top, z_bottom, &
                                       terrain_detected, zFaceProfile)
      call mpas_pool_get_config(domain % blocklist % configs, 'config_electrostatic_terrain_mode', terrain_mode)

      zgrid_mode = .false.
      if (terrain_detected) then
         select case (trim(terrain_mode))
         case ('zgrid')
            zgrid_mode = .true.
         case default
            call mpas_log_write('Electrostatics requires a flat (z-level) mesh: layer ' // &
                                'thickness varies horizontally, indicating a terrain-following ' // &
                                'coordinate. Set config_electrostatic_terrain_mode = ''zgrid'' ' // &
                                'to solve on a uniform altitude grid, or disable ' // &
                                'config_electrostatic_enable.', messageType=MPAS_LOG_CRIT)
         end select
      end if

      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         allocate(new_state)
         nullify(new_state % next)

         call init_block_state(block_ptr, new_state, domain % sphere_radius, source_x0, source_y0, &
                               mms_x, mms_y, z_top, z_bottom, zFaceProfile, zgrid_mode)

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
         call prepare_block_potential_rhs(block_ptr, state_ptr)
         block_ptr => block_ptr % next
         state_ptr => state_ptr % next
      end do

      call solve_domain_potential(domain)
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
         if (allocated(state_ptr % x)) deallocate(state_ptr % x)
         if (allocated(state_ptr % zMid)) deallocate(state_ptr % zMid)
         if (allocated(state_ptr % zGround)) deallocate(state_ptr % zGround)
         if (allocated(state_ptr % zFaceP)) deallocate(state_ptr % zFaceP)
         if (allocated(state_ptr % zMidPC)) deallocate(state_ptr % zMidPC)
         if (allocated(state_ptr % airThickness)) deallocate(state_ptr % airThickness)
         if (allocated(state_ptr % wallWeight)) deallocate(state_ptr % wallWeight)
         if (allocated(state_ptr % work)) deallocate(state_ptr % work)
         if (allocated(state_ptr % kFirstActive)) deallocate(state_ptr % kFirstActive)
         if (allocated(state_ptr % active)) deallocate(state_ptr % active)
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
!  routine compute_domain_synthetic_source_center
!
!> \brief Compute a mesh-wide center for synthetic charge sources
!
!-----------------------------------------------------------------------
   subroutine compute_domain_synthetic_source_center(domain, source_x0, source_y0)

      type(domain_type), intent(inout) :: domain
      real(kind=RKIND), intent(out) :: source_x0
      real(kind=RKIND), intent(out) :: source_y0

      type(block_type), pointer :: block_ptr
      type(mpas_pool_type), pointer :: mesh_pool
      integer, pointer :: nCellsSolve_ptr
      integer :: nCellsSolve
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell
      real(kind=RKIND) :: local_xmin, local_xmax, local_ymin, local_ymax
      real(kind=RKIND) :: global_xmin, global_xmax, global_ymin, global_ymax

      local_xmin = huge(1.0_RKIND)
      local_xmax = -huge(1.0_RKIND)
      local_ymin = huge(1.0_RKIND)
      local_ymax = -huge(1.0_RKIND)

      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
         call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve_ptr)
         call mpas_pool_get_array(mesh_pool, 'xCell', xCell)
         call mpas_pool_get_array(mesh_pool, 'yCell', yCell)

         nCellsSolve = nCellsSolve_ptr
         if (nCellsSolve > 0) then
            local_xmin = min(local_xmin, minval(xCell(1:nCellsSolve)))
            local_xmax = max(local_xmax, maxval(xCell(1:nCellsSolve)))
            local_ymin = min(local_ymin, minval(yCell(1:nCellsSolve)))
            local_ymax = max(local_ymax, maxval(yCell(1:nCellsSolve)))
         end if

         block_ptr => block_ptr % next
      end do

      call mpas_dmpar_min_real(domain % dminfo, local_xmin, global_xmin)
      call mpas_dmpar_max_real(domain % dminfo, local_xmax, global_xmax)
      call mpas_dmpar_min_real(domain % dminfo, local_ymin, global_ymin)
      call mpas_dmpar_max_real(domain % dminfo, local_ymax, global_ymax)

      if (global_xmax >= global_xmin) then
         source_x0 = 0.5_RKIND * (global_xmin + global_xmax)
      else
         source_x0 = 0.0_RKIND
      end if
      if (global_ymax >= global_ymin) then
         source_y0 = 0.5_RKIND * (global_ymin + global_ymax)
      else
         source_y0 = 0.0_RKIND
      end if

   end subroutine compute_domain_synthetic_source_center

!***********************************************************************
!
!  routine compute_domain_grid_extents
!
!> \brief Compute domain extents and globally detect terrain-following layers
!
!-----------------------------------------------------------------------
   subroutine compute_domain_grid_extents(domain, x_extent, y_extent, z_top, z_bottom, &
                                          terrain_detected, zFaceProfile)

      type(domain_type), intent(inout) :: domain
      real(kind=RKIND), intent(out) :: x_extent
      real(kind=RKIND), intent(out) :: y_extent
      real(kind=RKIND), intent(out) :: z_top
      real(kind=RKIND), intent(out) :: z_bottom
      logical, intent(out) :: terrain_detected
      real(kind=RKIND), dimension(:), allocatable, intent(out) :: zFaceProfile

      type(block_type), pointer :: block_ptr
      type(mpas_pool_type), pointer :: mesh_pool
      integer, pointer :: nCellsSolve_ptr, nVertLevels_ptr
      integer :: nCellsSolve, nVertLevels, k, local_flag, global_flag, nProfile
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell
      real(kind=RKIND), dimension(:,:), pointer :: zgrid
      real(kind=RKIND) :: lxmin, lxmax, lymin, lymax, lzmax, lzgmin
      real(kind=RKIND) :: gxmin, gxmax, gymin, gymax
      real(kind=RKIND) :: mean_dz, max_dev
      real(kind=RKIND), allocatable :: profile_local(:), profile_global(:)
      real(kind=RKIND), parameter :: flat_rtol = 1.0e-3_RKIND

      lxmin = huge(1.0_RKIND)
      lxmax = -huge(1.0_RKIND)
      lymin = huge(1.0_RKIND)
      lymax = -huge(1.0_RKIND)
      lzmax = -huge(1.0_RKIND)
      lzgmin = huge(1.0_RKIND)
      local_flag = 0

      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
         call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve_ptr)
         call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels_ptr)
         call mpas_pool_get_array(mesh_pool, 'xCell', xCell)
         call mpas_pool_get_array(mesh_pool, 'yCell', yCell)
         call mpas_pool_get_array(mesh_pool, 'zgrid', zgrid)

         nCellsSolve = nCellsSolve_ptr
         nVertLevels = nVertLevels_ptr

         if (.not. allocated(profile_local)) then
            allocate(profile_local(nVertLevels+2))   ! K+1 face sums + owned-cell count
            profile_local = 0.0_RKIND
         end if

         if (nCellsSolve > 0) then
            lxmin = min(lxmin, minval(xCell(1:nCellsSolve)))
            lxmax = max(lxmax, maxval(xCell(1:nCellsSolve)))
            lymin = min(lymin, minval(yCell(1:nCellsSolve)))
            lymax = max(lymax, maxval(yCell(1:nCellsSolve)))
            lzmax = max(lzmax, maxval(zgrid(nVertLevels+1,1:nCellsSolve)))
            lzgmin = min(lzgmin, minval(zgrid(1,1:nCellsSolve)))

            do k = 1, nVertLevels
               mean_dz = sum(zgrid(k+1,1:nCellsSolve) - zgrid(k,1:nCellsSolve)) / &
                         real(nCellsSolve, kind=RKIND)
               max_dev = maxval(abs((zgrid(k+1,1:nCellsSolve) - zgrid(k,1:nCellsSolve)) - mean_dz))
               if (mean_dz > 0.0_RKIND .and. max_dev > flat_rtol * mean_dz) local_flag = 1
            end do

            do k = 1, nVertLevels + 1
               profile_local(k) = profile_local(k) + sum(zgrid(k,1:nCellsSolve))
            end do
            profile_local(nVertLevels+2) = profile_local(nVertLevels+2) + &
                                           real(nCellsSolve, kind=RKIND)
         end if

         block_ptr => block_ptr % next
      end do

      call mpas_dmpar_min_real(domain % dminfo, lxmin, gxmin)
      call mpas_dmpar_max_real(domain % dminfo, lxmax, gxmax)
      call mpas_dmpar_min_real(domain % dminfo, lymin, gymin)
      call mpas_dmpar_max_real(domain % dminfo, lymax, gymax)
      call mpas_dmpar_max_real(domain % dminfo, lzmax, z_top)
      call mpas_dmpar_min_real(domain % dminfo, lzgmin, z_bottom)
      call mpas_dmpar_max_int(domain % dminfo, local_flag, global_flag)
      terrain_detected = global_flag > 0

      x_extent = gxmax - gxmin
      y_extent = gymax - gymin
      if (x_extent <= 0.0_RKIND) x_extent = 1.0_RKIND
      if (y_extent <= 0.0_RKIND) y_extent = 1.0_RKIND
      if (z_top <= z_bottom) z_top = z_bottom + 1.0_RKIND

      ! Global per-level mean interface profile (flat meshes adopt these faces
      ! natively). A single array allreduce keeps every rank bit-identical --
      ! unlike the retired per-rank layer-thickness means, which could disagree
      ! across ranks at roundoff level. Assumes >= 1 block per rank (standard
      ! MPAS decomposition); nVertLevels is a registry dimension, identical
      ! on all ranks.
      nProfile = size(profile_local)
      allocate(profile_global(nProfile))
      call mpas_dmpar_sum_real_array(domain % dminfo, nProfile, profile_local, profile_global)
      allocate(zFaceProfile(nProfile-1))
      if (profile_global(nProfile) > 0.0_RKIND) then
         zFaceProfile(:) = profile_global(1:nProfile-1) / profile_global(nProfile)
      else
         do k = 1, nProfile - 1
            zFaceProfile(k) = z_bottom + real(k-1, kind=RKIND) * (z_top - z_bottom) / &
                              real(nProfile-2, kind=RKIND)
         end do
      end if
      deallocate(profile_local, profile_global)

   end subroutine compute_domain_grid_extents

!***********************************************************************
!
!  routine init_block_state
!
!> \brief Initialize one block's electrostatic operator and PCG workspace
!
!-----------------------------------------------------------------------
   subroutine init_block_state(block_ptr, state_ptr, sphere_radius, source_x0, source_y0, &
                               x_period_in, y_period_in, z_top_in, z_bottom_in, zFaceProfile, zgrid_mode)

      type(block_type), pointer, intent(inout) :: block_ptr
      type(electrostatic_block_state), pointer, intent(inout) :: state_ptr
      real(kind=RKIND), intent(in) :: sphere_radius
      real(kind=RKIND), intent(in) :: source_x0
      real(kind=RKIND), intent(in) :: source_y0
      real(kind=RKIND), intent(in) :: x_period_in
      real(kind=RKIND), intent(in) :: y_period_in
      real(kind=RKIND), intent(in) :: z_top_in
      real(kind=RKIND), intent(in) :: z_bottom_in
      real(kind=RKIND), dimension(:), intent(in) :: zFaceProfile
      logical, intent(in) :: zgrid_mode

      type(mpas_pool_type), pointer :: mesh_pool, es_pool
      integer, pointer :: nCells_ptr, nEdges_ptr, nVertLevels_ptr, maxEdges_ptr
      integer, pointer :: nCellsSolve_ptr
      integer, pointer :: max_iter
      integer :: nCells, nEdges, nVertLevels, maxEdges, preconditioner, k, nDiagCells
      real(kind=RKIND) :: dzP
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, cellsOnEdge
      real(kind=RKIND), dimension(:), pointer :: areaCell, dcEdge, dvEdge
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell
      real(kind=RKIND), dimension(:,:), pointer :: zgrid
      real(kind=RKIND), dimension(:,:), pointer :: phi, rho_charge, E_normal
      real(kind=RKIND), dimension(:,:), pointer :: poisson_zmid, poisson_air_thickness
      real(kind=RKIND), dimension(:,:,:), pointer :: E_vector
      real(kind=RKIND), pointer :: tol
      real(kind=RKIND), pointer :: mesh_x_period, mesh_y_period
      real(kind=RKIND), allocatable :: areaMetric(:), dcMetric(:), dvMetric(:)
      real(kind=RKIND), allocatable :: h_weight(:,:), v_up(:,:), v_lo(:,:)
      real(kind=RKIND), allocatable :: diag_extra(:,:), op_volume(:,:)
      logical, pointer :: is_periodic, on_a_sphere
      character(len=StrKIND), pointer :: pc_name

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
      call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)

      call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nEdges', nEdges_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'maxEdges', maxEdges_ptr)
      nCells = nCells_ptr
      state_ptr % nCellsSolve = nCellsSolve_ptr
      nEdges = nEdges_ptr
      nVertLevels = nVertLevels_ptr
      maxEdges = maxEdges_ptr

      call mpas_pool_get_array(mesh_pool, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh_pool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh_pool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh_pool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh_pool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh_pool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh_pool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh_pool, 'xCell', xCell)
      call mpas_pool_get_array(mesh_pool, 'yCell', yCell)
      call mpas_pool_get_array(mesh_pool, 'zgrid', zgrid)
      call mpas_pool_get_config(mesh_pool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(mesh_pool, 'is_periodic', is_periodic)
      call mpas_pool_get_config(mesh_pool, 'x_period', mesh_x_period)
      call mpas_pool_get_config(mesh_pool, 'y_period', mesh_y_period)

      call mpas_pool_get_config(block_ptr % configs, 'config_poisson_preconditioner', pc_name)
      call mpas_pool_get_config(block_ptr % configs, 'config_poisson_tol', tol)
      call mpas_pool_get_config(block_ptr % configs, 'config_poisson_max_iter', max_iter)

      allocate(state_ptr % rhs(nVertLevels,nCells))
      allocate(state_ptr % x(nVertLevels,nCells))
      allocate(state_ptr % zMid(nVertLevels,nCells))
      allocate(state_ptr % zGround(nCells))
      allocate(state_ptr % zFaceP(nVertLevels+1))
      allocate(state_ptr % zMidPC(nVertLevels,nCells))
      allocate(state_ptr % airThickness(nVertLevels,nCells))
      allocate(state_ptr % wallWeight(nVertLevels,nCells))
      allocate(state_ptr % work(nVertLevels,nCells))
      allocate(state_ptr % kFirstActive(nCells))
      allocate(state_ptr % active(nVertLevels,nCells))
      allocate(state_ptr % E_vertical(nVertLevels,nCells))
      allocate(state_ptr % E_x(nVertLevels,nCells))
      allocate(state_ptr % E_y(nVertLevels,nCells))
      allocate(state_ptr % E_z(nVertLevels,nCells))
      allocate(state_ptr % E_zonal(nVertLevels,nCells))
      allocate(state_ptr % E_meridional(nVertLevels,nCells))
      allocate(areaMetric(nCells))
      allocate(dcMetric(nEdges))
      allocate(dvMetric(nEdges))
      allocate(h_weight(nVertLevels,nEdges))
      allocate(v_up(nVertLevels,nCells), v_lo(nVertLevels,nCells))
      allocate(diag_extra(nVertLevels,nCells), op_volume(nVertLevels,nCells))

      state_ptr % zgrid_mode = zgrid_mode
      call fill_z_mid(nCells, nVertLevels, zgrid, state_ptr % zMid)
      state_ptr % zGround(1:nCells) = zgrid(1,1:nCells)

      if (zgrid_mode) then
         ! Terrain: uniform Poisson faces spanning [z_bottom, z_top] (the
         ! Phase-4 z-grid design; sources and outputs are remapped).
         dzP = (z_top_in - z_bottom_in) / real(nVertLevels, kind=RKIND)
         do k = 1, nVertLevels + 1
            state_ptr % zFaceP(k) = z_bottom_in + real(k-1, kind=RKIND) * dzP
         end do
      else
         ! Flat: adopt the native (horizontally uniform) MPAS interface
         ! profile so the operator, sources, and outputs share one grid.
         ! With zGround(i) = zgrid(1,i) = zFaceP(1) the builders give
         !   kFirstActive = 1, active = .true., airThickness(k) = dz(k),
         !   zMidPC(k) = 0.5*(zFace(k)+zFace(k+1)) = zMid(k),
         ! and the FV weights reduce to the legacy flat operator exactly
         ! (mpas_elliptic_operator.F flat wrapper):
         !   v_up(k) = eps*A/(zMidPC(k+1)-zMidPC(k)) = eps*A/(0.5*(dz(k)+dz(k+1)))
         !   v_lo(1) = eps*A/(zMidPC(1)-zGround)     = eps*A/(dz(1)/2)
         !   h(k,e)  = eps*(l_e/d_e)*dz(k),  volume(k) = A*dz(k),  no walls.
         ! Sources evaluated at zMid, rhs = rho*volume, and phi/E output at
         ! zMid are then mutually consistent on vertically stretched flat
         ! meshes -- no remap needed, and zMidPC == zMid by construction.
         state_ptr % zFaceP(1:nVertLevels+1) = zFaceProfile(1:nVertLevels+1)
      end if

      call electrostatic_build_poisson_grid(nCells, nVertLevels, state_ptr % zGround, &
                                            state_ptr % zFaceP, state_ptr % zMidPC, &
                                            state_ptr % active, state_ptr % kFirstActive, &
                                            state_ptr % airThickness)

      if (sphere_radius > 0.0_RKIND) state_ptr % sphere_radius = sphere_radius
      state_ptr % source_x0 = source_x0
      state_ptr % source_y0 = source_y0
      state_ptr % x_period = x_period_in
      state_ptr % y_period = y_period_in
      state_ptr % z_top = z_top_in
      state_ptr % horizontal_length_scale = &
         electrostatic_horizontal_metric_scale(nEdges, dcEdge, state_ptr % sphere_radius, on_a_sphere)
      areaMetric = areaCell * state_ptr % horizontal_length_scale * state_ptr % horizontal_length_scale
      dcMetric = dcEdge * state_ptr % horizontal_length_scale
      dvMetric = dvEdge * state_ptr % horizontal_length_scale
      if (is_periodic) then
         if (mesh_x_period > 0.0_RKIND) state_ptr % x_period = mesh_x_period
         if (mesh_y_period > 0.0_RKIND) state_ptr % y_period = mesh_y_period
      end if

      preconditioner = preconditioner_from_name(pc_name)
      call electrostatic_build_poisson_fv_weights(nCells, nEdges, nVertLevels, cellsOnEdge, &
                                                  areaMetric, dcMetric, dvMetric, epsilon0, &
                                                  state_ptr % zFaceP, state_ptr % zGround, &
                                                  state_ptr % z_top, state_ptr % zMidPC, &
                                                  state_ptr % active, state_ptr % kFirstActive, &
                                                  state_ptr % airThickness, h_weight, v_up, v_lo, &
                                                  diag_extra, op_volume)
      state_ptr % wallWeight = diag_extra
      call elliptic_operator_init_from_weights(state_ptr % op, nCells, nEdges, nVertLevels, maxEdges, &
                                               nEdgesOnCell, cellsOnCell, edgesOnCell, &
                                               h_weight, v_up, v_lo, diag_extra, op_volume, &
                                               state_ptr % active)
      call pcg_init(state_ptr % solver, nCells, nVertLevels, preconditioner, max(1, max_iter), tol)

      state_ptr % rhs = 0.0_RKIND
      state_ptr % x = 0.0_RKIND
      state_ptr % E_vertical = 0.0_RKIND
      state_ptr % E_x = 0.0_RKIND
      state_ptr % E_y = 0.0_RKIND
      state_ptr % E_z = 0.0_RKIND
      state_ptr % E_zonal = 0.0_RKIND
      state_ptr % E_meridional = 0.0_RKIND
      state_ptr % work = 0.0_RKIND

      call mpas_pool_get_array(es_pool, 'phi', phi)
      call mpas_pool_get_array(es_pool, 'rho_charge', rho_charge)
      call mpas_pool_get_array(es_pool, 'E_normal', E_normal)
      call mpas_pool_get_array(es_pool, 'E_vector', E_vector)
      call mpas_pool_get_array(es_pool, 'poisson_zmid', poisson_zmid)
      call mpas_pool_get_array(es_pool, 'poisson_air_thickness', poisson_air_thickness)
      phi = 0.0_RKIND
      rho_charge = 0.0_RKIND
      E_normal = 0.0_RKIND
      E_vector = 0.0_RKIND

      ! Poisson-grid diagnostics are time-invariant geometry (written once by
      ! the grid builder above); fill the registry fields here rather than
      ! every solve. The min clamp is load-bearing -- pool arrays carry the
      ! +1 garbage column.
      nDiagCells = min(size(poisson_zmid, 2), size(state_ptr % zMidPC, 2))
      poisson_zmid = 0.0_RKIND
      poisson_air_thickness = 0.0_RKIND
      poisson_zmid(:,1:nDiagCells) = state_ptr % zMidPC(:,1:nDiagCells)
      poisson_air_thickness(:,1:nDiagCells) = state_ptr % airThickness(:,1:nDiagCells)

      call set_solver_diagnostics(es_pool, 0, 0.0_RKIND, 0.0_RKIND)

      deallocate(areaMetric)
      deallocate(dcMetric)
      deallocate(dvMetric)
      deallocate(h_weight, v_up, v_lo, diag_extra, op_volume)

   end subroutine init_block_state

!***********************************************************************
!
!  routine prepare_block_potential_rhs
!
!> \brief Fill source terms and RHS for one block
!
!-----------------------------------------------------------------------
   subroutine prepare_block_potential_rhs(block_ptr, state_ptr)

      type(block_type), pointer, intent(inout) :: block_ptr
      type(electrostatic_block_state), pointer, intent(inout) :: state_ptr

      type(mpas_pool_type), pointer :: mesh_pool, es_pool, state_pool
      integer, pointer :: nCells_ptr, nVertLevels_ptr
      integer, pointer :: num_scalars_ptr, index_qc_ptr, index_qr_ptr, index_qi_ptr, index_qg_ptr
      integer :: nCells, nVertLevels, stub_negative_index, stub_positive_index
      real(kind=RKIND), dimension(:), pointer :: latCell, lonCell, xCell, yCell
      real(kind=RKIND), dimension(:,:), pointer :: zgrid
      real(kind=RKIND), dimension(:,:), pointer :: phi, rho_charge
      real(kind=RKIND), dimension(:,:), pointer :: w
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars
      real(kind=RKIND), pointer :: phi_ground
      real(kind=RKIND), pointer :: source_charge, source_sigma, dipole_separation
      real(kind=RKIND), pointer :: hill_height
      real(kind=RKIND), pointer :: stub_alpha, stub_beta
      character(len=StrKIND), pointer :: source_name

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh_pool)
      call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', state_pool)

      call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells_ptr)
      call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels_ptr)
      call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars_ptr)
      call mpas_pool_get_dimension(state_pool, 'index_qc', index_qc_ptr)
      call mpas_pool_get_dimension(state_pool, 'index_qr', index_qr_ptr)
      call mpas_pool_get_dimension(state_pool, 'index_qi', index_qi_ptr)
      call mpas_pool_get_dimension(state_pool, 'index_qg', index_qg_ptr)
      nCells = nCells_ptr
      nVertLevels = nVertLevels_ptr

      call mpas_pool_get_array(mesh_pool, 'xCell', xCell)
      call mpas_pool_get_array(mesh_pool, 'yCell', yCell)
      call mpas_pool_get_array(mesh_pool, 'latCell', latCell)
      call mpas_pool_get_array(mesh_pool, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh_pool, 'zgrid', zgrid)
      call mpas_pool_get_array(state_pool, 'w', w, 1)
      call mpas_pool_get_array(state_pool, 'scalars', scalars, 1)
      call mpas_pool_get_array(es_pool, 'phi', phi)
      call mpas_pool_get_array(es_pool, 'rho_charge', rho_charge)

      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_source', source_name)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_bc_ground', phi_ground)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_source_charge', source_charge)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_source_sigma', source_sigma)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_dipole_separation', dipole_separation)
      call mpas_pool_get_config(block_ptr % configs, 'config_electrostatic_hill_height', hill_height)
      call mpas_pool_get_config(block_ptr % configs, 'config_stub_alpha', stub_alpha)
      call mpas_pool_get_config(block_ptr % configs, 'config_stub_beta', stub_beta)

      call fill_z_mid(nCells, nVertLevels, zgrid, state_ptr % zMid)

      if (state_ptr % zgrid_mode) then
         select case (trim(source_name))
         case ('mms_terrain', 'stub', 'external')
            continue
         case default
            call mpas_log_write('Electrostatics: source "'//trim(source_name)// &
               '" is not supported with config_electrostatic_terrain_mode=zgrid; ' // &
               'use mms_terrain, stub, or external.', messageType=MPAS_LOG_CRIT)
         end select
      end if

      select case (trim(source_name))
      case ('external')
         continue
      case ('mms_terrain')
         if (.not. state_ptr % zgrid_mode) then
            call mpas_log_write('Electrostatics: source mms_terrain requires a terrain ' // &
               'mesh with config_electrostatic_terrain_mode = ''zgrid''.', &
               messageType=MPAS_LOG_CRIT)
         else
            call electrostatic_fill_mms_terrain_source(nCells, nVertLevels, xCell, yCell, &
                                                       state_ptr % zMidPC, epsilon0, &
                                                       state_ptr % x_period, state_ptr % y_period, &
                                                       state_ptr % z_top, hill_height, &
                                                       state_ptr % work, active=state_ptr % active)
            call electrostatic_fill_mms_terrain_source(nCells, nVertLevels, xCell, yCell, &
                                                       state_ptr % zMid, epsilon0, &
                                                       state_ptr % x_period, state_ptr % y_period, &
                                                       state_ptr % z_top, hill_height, rho_charge)
         end if
      case ('mms_cart_horizontal')
         call electrostatic_fill_mms_cart_horizontal_source(nCells, nVertLevels, xCell, yCell, &
                                                            state_ptr % zMid, epsilon0, &
                                                            state_ptr % x_period, state_ptr % y_period, &
                                                            state_ptr % z_top, &
                                                            state_ptr % op % v_weight_upper, &
                                                            state_ptr % op % v_weight_lower, &
                                                            state_ptr % op % volume, rho_charge)
      case ('mms_sphere')
         call electrostatic_fill_mms_sphere_source(nCells, nVertLevels, latCell, lonCell, &
                                                   state_ptr % zMid, epsilon0, &
                                                   state_ptr % sphere_radius, state_ptr % z_top, &
                                                   rho_charge)
      case ('mms_sphere_horizontal')
         call electrostatic_fill_mms_sphere_horizontal_source(nCells, nVertLevels, latCell, lonCell, &
                                                              state_ptr % zMid, epsilon0, &
                                                              state_ptr % sphere_radius, state_ptr % z_top, &
                                                              state_ptr % op % v_weight_upper, &
                                                              state_ptr % op % v_weight_lower, &
                                                              state_ptr % op % volume, rho_charge)
      case ('dipole')
         call electrostatic_fill_dipole_source(nCells, nVertLevels, xCell, yCell, &
                                               state_ptr % zMid, state_ptr % source_x0, &
                                               state_ptr % source_y0, rho_charge)
      case ('tripole')
         call electrostatic_fill_tripole_source(nCells, nVertLevels, xCell, yCell, &
                                                state_ptr % zMid, state_ptr % source_x0, &
                                                state_ptr % source_y0, rho_charge)
      case ('gaussian_monopole')
         call electrostatic_fill_gaussian_monopole_source(nCells, nVertLevels, xCell, yCell, &
                                                          state_ptr % zMid, state_ptr % source_x0, &
                                                          state_ptr % source_y0, 0.5_RKIND * state_ptr % z_top, &
                                                          source_charge, source_sigma, rho_charge)
      case ('gaussian_dipole_y')
         call electrostatic_fill_gaussian_dipole_y_source(nCells, nVertLevels, xCell, yCell, &
                                                          state_ptr % zMid, state_ptr % source_x0, &
                                                          state_ptr % source_y0, 0.5_RKIND * state_ptr % z_top, &
                                                          source_charge, source_sigma, dipole_separation, &
                                                          rho_charge)
      case ('gaussian_dipole_z')
         call electrostatic_fill_gaussian_dipole_z_source(nCells, nVertLevels, xCell, yCell, &
                                                          state_ptr % zMid, state_ptr % source_x0, &
                                                          state_ptr % source_y0, 0.5_RKIND * state_ptr % z_top, &
                                                          source_charge, source_sigma, dipole_separation, &
                                                          rho_charge)
      case ('stub')
         call electrostatic_select_stub_indices(num_scalars_ptr, index_qi_ptr, index_qg_ptr, &
                                                index_qc_ptr, index_qr_ptr, &
                                                stub_negative_index, stub_positive_index)
         call electrostatic_fill_stub_source(nCells, nVertLevels, num_scalars_ptr, w, scalars, &
                                             stub_negative_index, stub_positive_index, stub_alpha, stub_beta, &
                                             rho_charge)
      case default
         call electrostatic_fill_mms_cart_source(nCells, nVertLevels, xCell, yCell, state_ptr % zMid, &
                                                 epsilon0, state_ptr % x_period, state_ptr % y_period, &
                                                 state_ptr % z_top, rho_charge)
      end select

      if (state_ptr % zgrid_mode) then
         select case (trim(source_name))
         case ('stub', 'external')
            call electrostatic_remap_rho_conservative(nCells, nVertLevels, zgrid, &
                     state_ptr % zFaceP, state_ptr % zGround, state_ptr % kFirstActive, &
                     state_ptr % active, state_ptr % airThickness, rho_charge, state_ptr % work)
         case default
            continue
         end select
         ! state%work holds the Poisson-grid charge density from the fill or
         ! remap above; it is consumed into rhs here, before the solve, and is
         ! reused as back-remap scratch only after the solve (postprocess).
         state_ptr % rhs = state_ptr % work * state_ptr % op % volume
      else
         state_ptr % rhs = rho_charge * state_ptr % op % volume
      end if

      if (abs(phi_ground) > 0.0_RKIND) then
         call electrostatic_apply_ground_rhs(nCells, nVertLevels, phi_ground, &
                                             state_ptr % op % v_weight_lower, state_ptr % rhs, &
                                             kFirstActive=state_ptr % kFirstActive, &
                                             wall_weight=state_ptr % wallWeight)
      end if

      if (state_ptr % zgrid_mode) then
         where (.not. state_ptr % active) state_ptr % x = 0.0_RKIND
      else
         state_ptr % x = phi
      end if

   end subroutine prepare_block_potential_rhs

!***********************************************************************
!
!  routine solve_domain_potential
!
!> \brief Solve A phi = rhs over the full decomposed domain
!
!-----------------------------------------------------------------------
   subroutine solve_domain_potential(domain)

      type(domain_type), intent(inout) :: domain

      type(electrostatic_block_state), pointer :: state_ptr
      integer :: iter, iter_count, max_iter, status
      real(kind=RKIND) :: alpha, beta, denom, rhs_norm, r0, rfinal, rz_new, rz_old, tol

      if (.not. associated(block_states)) return

      max_iter = max(1, block_states % solver % max_iter)
      tol = block_states % solver % tol
      iter_count = 0
      status = 1

      call matvec_domain_vector(domain, 'x')

      state_ptr => block_states
      do while (associated(state_ptr))
         state_ptr % solver % r = 0.0_RKIND
         state_ptr % solver % r(:,1:state_ptr % nCellsSolve) = &
            state_ptr % rhs(:,1:state_ptr % nCellsSolve) - &
            state_ptr % solver % Ap(:,1:state_ptr % nCellsSolve)
         state_ptr => state_ptr % next
      end do

      rhs_norm = sqrt(domain_dot('rhs', 'rhs', domain))
      r0 = domain_residual_norm('r', rhs_norm, domain)
      rfinal = r0

      if (rfinal <= tol) then
         status = 0
         call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
         return
      end if

      if (.not. apply_domain_preconditioner(status)) then
         call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
         return
      end if

      rz_old = domain_dot('r', 'z', domain)
      if (rz_old <= 0.0_RKIND) then
         status = 2
         call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
         return
      end if

      state_ptr => block_states
      do while (associated(state_ptr))
         state_ptr % solver % p = 0.0_RKIND
         state_ptr % solver % p(:,1:state_ptr % nCellsSolve) = &
            state_ptr % solver % z(:,1:state_ptr % nCellsSolve)
         state_ptr => state_ptr % next
      end do

      do iter = 1, max_iter
         call matvec_domain_vector(domain, 'p')
         denom = domain_dot('p', 'Ap', domain)
         if (denom <= 0.0_RKIND) then
            iter_count = iter - 1
            status = 2
            call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
            return
         end if

         alpha = rz_old / denom
         state_ptr => block_states
         do while (associated(state_ptr))
            state_ptr % x(:,1:state_ptr % nCellsSolve) = &
               state_ptr % x(:,1:state_ptr % nCellsSolve) + &
               alpha * state_ptr % solver % p(:,1:state_ptr % nCellsSolve)
            state_ptr % solver % r(:,1:state_ptr % nCellsSolve) = &
               state_ptr % solver % r(:,1:state_ptr % nCellsSolve) - &
               alpha * state_ptr % solver % Ap(:,1:state_ptr % nCellsSolve)
            state_ptr => state_ptr % next
         end do

         rfinal = domain_residual_norm('r', rhs_norm, domain)
         iter_count = iter
         if (rfinal <= tol) then
            status = 0
            call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
            return
         end if

         if (.not. apply_domain_preconditioner(status)) then
            call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
            return
         end if

         rz_new = domain_dot('r', 'z', domain)
         if (rz_new <= 0.0_RKIND) then
            status = 2
            call finish_domain_pcg(domain, iter_count, r0, rfinal, status)
            return
         end if

         beta = rz_new / rz_old
         state_ptr => block_states
         do while (associated(state_ptr))
            state_ptr % solver % p(:,1:state_ptr % nCellsSolve) = &
               state_ptr % solver % z(:,1:state_ptr % nCellsSolve) + &
               beta * state_ptr % solver % p(:,1:state_ptr % nCellsSolve)
            state_ptr => state_ptr % next
         end do
         rz_old = rz_new
      end do

      status = 1
      call finish_domain_pcg(domain, iter_count, r0, rfinal, status)

   end subroutine solve_domain_potential

!***********************************************************************
!
!  routine matvec_domain_vector
!
!> \brief Apply the operator to a decomposed vector with halo exchange
!
!-----------------------------------------------------------------------
   subroutine matvec_domain_vector(domain, vector_name)

      type(domain_type), intent(inout) :: domain
      character(len=*), intent(in) :: vector_name

      type(block_type), pointer :: block_ptr
      type(electrostatic_block_state), pointer :: state_ptr
      type(mpas_pool_type), pointer :: es_pool
      real(kind=RKIND), dimension(:,:), pointer :: phi

      block_ptr => domain % blocklist
      state_ptr => block_states
      do while (associated(block_ptr) .and. associated(state_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)
         call mpas_pool_get_array(es_pool, 'phi', phi)
         phi = 0.0_RKIND
         select case (trim(vector_name))
         case ('x')
            phi(:,1:state_ptr % nCellsSolve) = state_ptr % x(:,1:state_ptr % nCellsSolve)
         case ('p')
            phi(:,1:state_ptr % nCellsSolve) = state_ptr % solver % p(:,1:state_ptr % nCellsSolve)
         end select
         block_ptr => block_ptr % next
         state_ptr => state_ptr % next
      end do

      call exchange_halo_group(domain, 'electrostatic:phi')

      block_ptr => domain % blocklist
      state_ptr => block_states
      do while (associated(block_ptr) .and. associated(state_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)
         call mpas_pool_get_array(es_pool, 'phi', phi)
         call state_ptr % op % matvec(phi, state_ptr % solver % Ap)
         block_ptr => block_ptr % next
         state_ptr => state_ptr % next
      end do

   end subroutine matvec_domain_vector

!***********************************************************************
!
!  routine apply_domain_preconditioner
!
!> \brief Apply the configured block-Jacobi preconditioner to owned cells
!
!-----------------------------------------------------------------------
   logical function apply_domain_preconditioner(status)

      integer, intent(out) :: status

      type(electrostatic_block_state), pointer :: state_ptr

      apply_domain_preconditioner = .false.

      state_ptr => block_states
      do while (associated(state_ptr))
         state_ptr % solver % z = 0.0_RKIND
         select case (state_ptr % solver % preconditioner)
         case (PC_NONE)
            state_ptr % solver % z(:,1:state_ptr % nCellsSolve) = &
               state_ptr % solver % r(:,1:state_ptr % nCellsSolve)
         case (PC_JACOBI)
            if (any(state_ptr % op % diag(:,1:state_ptr % nCellsSolve) <= 0.0_RKIND)) then
               status = 3
               return
            end if
            state_ptr % solver % z(:,1:state_ptr % nCellsSolve) = &
               state_ptr % solver % r(:,1:state_ptr % nCellsSolve) / &
               state_ptr % op % diag(:,1:state_ptr % nCellsSolve)
         case default
            status = 3
            return
         end select
         state_ptr => state_ptr % next
      end do

      status = 1
      apply_domain_preconditioner = .true.

   end function apply_domain_preconditioner

!***********************************************************************
!
!  routine finish_domain_pcg
!
!> \brief Copy the converged PCG solution to phi and record diagnostics
!
!-----------------------------------------------------------------------
   subroutine finish_domain_pcg(domain, iter_count, r0, rfinal, status)

      type(domain_type), intent(inout) :: domain
      integer, intent(in) :: iter_count
      real(kind=RKIND), intent(in) :: r0
      real(kind=RKIND), intent(in) :: rfinal
      integer, intent(in) :: status

      type(block_type), pointer :: block_ptr
      type(electrostatic_block_state), pointer :: state_ptr
      type(mpas_pool_type), pointer :: es_pool
      real(kind=RKIND), dimension(:,:), pointer :: phi

      block_ptr => domain % blocklist
      state_ptr => block_states
      do while (associated(block_ptr) .and. associated(state_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'electrostatic', es_pool)
         call mpas_pool_get_array(es_pool, 'phi', phi)
         phi = 0.0_RKIND
         phi(:,1:state_ptr % nCellsSolve) = state_ptr % x(:,1:state_ptr % nCellsSolve)
         call set_solver_diagnostics(es_pool, iter_count, r0, rfinal)
         block_ptr => block_ptr % next
         state_ptr => state_ptr % next
      end do

      if (status /= 0) then
         call mpas_log_write('Electrostatic PCG returned status $i after $i iterations; final residual = $r', &
                             messageType=MPAS_LOG_WARN, intArgs=(/status, iter_count/), realArgs=(/rfinal/))
      end if

   end subroutine finish_domain_pcg

!***********************************************************************
!
!  routine domain_residual_norm
!
!> \brief Return a global relative residual norm
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function domain_residual_norm(vector_name, rhs_norm, domain)

      character(len=*), intent(in) :: vector_name
      real(kind=RKIND), intent(in) :: rhs_norm
      type(domain_type), intent(in) :: domain

      real(kind=RKIND) :: r_norm

      r_norm = sqrt(domain_dot(vector_name, vector_name, domain))
      if (rhs_norm > 0.0_RKIND) then
         domain_residual_norm = r_norm / rhs_norm
      else
         domain_residual_norm = r_norm
      end if

   end function domain_residual_norm

!***********************************************************************
!
!  routine domain_dot
!
!> \brief Return a global dot product over owned cells only
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function domain_dot(a_name, b_name, domain)

      character(len=*), intent(in) :: a_name
      character(len=*), intent(in) :: b_name
      type(domain_type), intent(in) :: domain

      type(electrostatic_block_state), pointer :: state_ptr
      real(kind=RKIND) :: local_dot, global_dot

      local_dot = 0.0_RKIND
      state_ptr => block_states
      do while (associated(state_ptr))
         local_dot = local_dot + block_dot(state_ptr, a_name, b_name)
         state_ptr => state_ptr % next
      end do

      call mpas_dmpar_sum_real(domain % dminfo, local_dot, global_dot)
      domain_dot = global_dot

   end function domain_dot

!***********************************************************************
!
!  routine block_dot
!
!> \brief Return a block-local dot product over owned cells only
!
!-----------------------------------------------------------------------
   real(kind=RKIND) function block_dot(state_ptr, a_name, b_name)

      type(electrostatic_block_state), pointer, intent(in) :: state_ptr
      character(len=*), intent(in) :: a_name
      character(len=*), intent(in) :: b_name

      select case (trim(a_name)//':'//trim(b_name))
      case ('rhs:rhs')
         block_dot = sum(state_ptr % rhs(:,1:state_ptr % nCellsSolve) * &
                         state_ptr % rhs(:,1:state_ptr % nCellsSolve))
      case ('r:r')
         block_dot = sum(state_ptr % solver % r(:,1:state_ptr % nCellsSolve) * &
                         state_ptr % solver % r(:,1:state_ptr % nCellsSolve))
      case ('r:z')
         block_dot = sum(state_ptr % solver % r(:,1:state_ptr % nCellsSolve) * &
                         state_ptr % solver % z(:,1:state_ptr % nCellsSolve))
      case ('p:Ap')
         block_dot = sum(state_ptr % solver % p(:,1:state_ptr % nCellsSolve) * &
                         state_ptr % solver % Ap(:,1:state_ptr % nCellsSolve))
      case default
         block_dot = 0.0_RKIND
      end select

   end function block_dot

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
      integer :: iCell, k, nCells, nEdges, nPhiCells, nVertLevels
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

      call electrostatic_compute_edge_normal_E(nCells, nEdges, nVertLevels, cellsOnEdge, dcEdge, &
                                               state_ptr % horizontal_length_scale, phi_ground, &
                                               state_ptr % active, phi, E_normal)

      call mpas_reconstruct(mesh_pool, E_normal, state_ptr % E_x, state_ptr % E_y, state_ptr % E_z, &
                            state_ptr % E_zonal, state_ptr % E_meridional, includeHalos=.false.)

      call fill_z_mid(nCells, nVertLevels, zgrid, state_ptr % zMid)
      if (state_ptr % zgrid_mode) then
         call electrostatic_compute_vertical_E(nCells, nVertLevels, phi, state_ptr % zMidPC, &
                                               phi_ground, state_ptr % z_top, state_ptr % zGround, &
                                               state_ptr % kFirstActive, state_ptr % E_vertical)

         call backremap_field_to_mpas(state_ptr, nCells, nVertLevels, state_ptr % E_x)
         call backremap_field_to_mpas(state_ptr, nCells, nVertLevels, state_ptr % E_y)
         call backremap_field_to_mpas(state_ptr, nCells, nVertLevels, state_ptr % E_z)
         call backremap_field_to_mpas(state_ptr, nCells, nVertLevels, state_ptr % E_vertical)
      else
         call electrostatic_compute_vertical_E(nCells, nVertLevels, phi, state_ptr % zMid, &
                                               phi_ground, state_ptr % z_top, state_ptr % zGround, &
                                               state_ptr % kFirstActive, state_ptr % E_vertical)
      end if

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

      if (state_ptr % zgrid_mode) then
         ! Registry phi is (nVertLevels, nCells+1) with a garbage column;
         ! passing it straight to the (nVertLevels,nCells) intent(in) source
         ! dummy sequence-associates the first nCells columns -- the same
         ! values the former reallocating phiP = phi copy produced.
         call electrostatic_interp_poisson_to_mpas(nCells, nVertLevels, state_ptr % zMidPC, &
                                                   state_ptr % zMid, state_ptr % zGround, &
                                                   state_ptr % kFirstActive, &
                                                   phi, state_ptr % work, lower_value=phi_ground)
         nPhiCells = min(size(phi, 2), size(state_ptr % work, 2))
         phi = 0.0_RKIND
         phi(:,1:nPhiCells) = state_ptr % work(:,1:nPhiCells)
      end if

   end subroutine postprocess_block_E

!***********************************************************************
!
!  routine backremap_field_to_mpas
!
!> \brief Remap one solver-grid field to MPAS midpoints in place
!
!-----------------------------------------------------------------------
   subroutine backremap_field_to_mpas(state_ptr, nCells, nVertLevels, field)

      type(electrostatic_block_state), pointer, intent(in) :: state_ptr
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      real(kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: field

      call electrostatic_interp_poisson_to_mpas(nCells, nVertLevels, state_ptr % zMidPC, &
                                                state_ptr % zMid, state_ptr % zGround, &
                                                state_ptr % kFirstActive, &
                                                field, state_ptr % work)
      field = state_ptr % work

   end subroutine backremap_field_to_mpas

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
      case ('jacobi')
         preconditioner_from_name = PC_JACOBI
      case default
         call mpas_log_write('Electrostatics: unrecognized config_poisson_preconditioner ' // &
                             '"'//trim(pc_name)//'"; falling back to Jacobi. Valid options are ' // &
                             '"none" and "jacobi" (block-SGS is not implemented in Phase 1).', &
                             messageType=MPAS_LOG_WARN)
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
