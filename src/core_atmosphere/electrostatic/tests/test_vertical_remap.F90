program test_vertical_remap

   use mpas_kind_types, only : RKIND
   use mpas_electrostatic_vertical_remap, only : electrostatic_interp_mpas_to_poisson, &
                                                electrostatic_interp_w_to_poisson_mid

   implicit none

   call test_linear_profile_exact()
   call test_inactive_poisson_levels_zero()
   call test_w_interfaces_interpolate_to_poisson_midpoints()
   print *, "PASS: vertical remap tests"

contains

   subroutine test_linear_profile_exact()
      integer, parameter :: nCells = 1, nVertLevels = 4
      real(kind=RKIND) :: zMpas(nVertLevels,nCells), zPoisson(nVertLevels,nCells)
      real(kind=RKIND) :: source(nVertLevels,nCells), target(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      zMpas(:,1) = [100.0_RKIND, 300.0_RKIND, 700.0_RKIND, 900.0_RKIND]
      source(:,1) = 2.0_RKIND * zMpas(:,1) + 5.0_RKIND
      zPoisson(:,1) = [100.0_RKIND, 200.0_RKIND, 500.0_RKIND, 900.0_RKIND]
      active(:,1) = .true.

      call electrostatic_interp_mpas_to_poisson(nCells, nVertLevels, zMpas, zPoisson, active, source, target)

      if (maxval(abs(target(:,1) - (2.0_RKIND*zPoisson(:,1) + 5.0_RKIND))) > 1.0e-12_RKIND) &
         stop "FAIL: linear interpolation not exact"
   end subroutine test_linear_profile_exact

   subroutine test_inactive_poisson_levels_zero()
      integer, parameter :: nCells = 1, nVertLevels = 3
      real(kind=RKIND) :: zMpas(nVertLevels,nCells), zPoisson(nVertLevels,nCells)
      real(kind=RKIND) :: source(nVertLevels,nCells), target(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      zMpas(:,1) = [100.0_RKIND, 300.0_RKIND, 500.0_RKIND]
      source(:,1) = [10.0_RKIND, 30.0_RKIND, 50.0_RKIND]
      zPoisson(:,1) = [50.0_RKIND, 100.0_RKIND, 200.0_RKIND]
      active(:,1) = [.false., .true., .true.]

      call electrostatic_interp_mpas_to_poisson(nCells, nVertLevels, zMpas, zPoisson, active, source, target)

      if (target(1,1) /= 0.0_RKIND) stop "FAIL: inactive level not zeroed"
      if (abs(target(2,1) - 10.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: lower clamp mismatch"
      if (abs(target(3,1) - 20.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: interior interpolation mismatch"
   end subroutine test_inactive_poisson_levels_zero

   subroutine test_w_interfaces_interpolate_to_poisson_midpoints()
      integer, parameter :: nCells = 1, nVertLevels = 2
      real(kind=RKIND) :: zgrid(nVertLevels+1,nCells), zPoisson(nVertLevels,nCells)
      real(kind=RKIND) :: w(nVertLevels+1,nCells), wPoisson(nVertLevels,nCells)
      logical :: active(nVertLevels,nCells)

      zgrid(:,1) = [0.0_RKIND, 100.0_RKIND, 300.0_RKIND]
      w(:,1) = [0.0_RKIND, 10.0_RKIND, 30.0_RKIND]
      zPoisson(:,1) = [50.0_RKIND, 200.0_RKIND]
      active(:,1) = .true.

      call electrostatic_interp_w_to_poisson_mid(nCells, nVertLevels, zgrid, zPoisson, active, w, wPoisson)

      if (abs(wPoisson(1,1) - 5.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: w interpolation level 1"
      if (abs(wPoisson(2,1) - 20.0_RKIND) > 1.0e-12_RKIND) stop "FAIL: w interpolation level 2"
   end subroutine test_w_interfaces_interpolate_to_poisson_midpoints

end program test_vertical_remap
