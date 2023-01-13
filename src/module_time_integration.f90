module time_integrator
  use field, only: field_type
  implicit none

  type, abstract :: time_integrator_type
     real :: starttime, endtime, dt
   contains
     procedure(integrate_proc), public, &
          & deferred :: integrate
  end type time_integrator_type

  type, extends(time_integrator_type) :: euler_integrator_type
     real :: alpha = 1.
   contains
     procedure :: integrate => integrate_euler
  end type euler_integrator_type

  type, extends(time_integrator_type) :: AB2_integrator_type
   contains
     procedure :: integrate => integrate_AB2
  end type AB2_integrator_type

  type, extends(time_integrator_type) :: RK3_integrator_type
   contains
     procedure :: integrate => integrate_RK3
  end type RK3_integrator_type

  abstract interface
     subroutine integrate_proc(self, afield)
       import time_integrator_type, field_type
       class(time_integrator_type), intent(in) :: self
       class(field_type), intent(inout) :: afield
     end subroutine integrate_proc
  end interface

contains

  subroutine integrate_euler(self, afield)
    class(euler_integrator_type), intent(in) :: self
    class(field_type), allocatable, intent(inout) :: afield
    integer :: nt
    integer :: i

    nt = floor((self%endtime - self%starttime) / self%dt)
    do i = 1, nt
       afield = euler_timestep(afield, self%dt)
    end do
  end subroutine integrate_euler

  subroutine integrate_AB2(self, afield)
    class(AB2_integrator_type), intent(in) :: self
    class(field_type), allocatable, intent(inout) :: afield
    class(field_type), allocatable :: f1, f2
    integer :: nt
    integer :: i

    nt = floor((self%endtime - self%starttime) / self%dt)
    f1 = afield
    f2 = euler_timestep(afield, self%dt)
    do i = 2, nt
       call AB2_timestep(f1, f2, self%dt)
    end do
    afield = f2
  end subroutine integrate_AB2

  subroutine integrate_RK3(self, afield)
    class(RK3_integrator_type), intent(in) :: self
    class(field_type), allocatable, intent(inout) :: afield
    class(field_type), allocatable :: afield2
    integer :: nt, i

    nt = floor((self%endtime - self%starttime) / self%dt)
    do i = 1, nt
       ! First fractional step is Euler like
       afield2 = afield%rhs() * (8. / 15.) * self%dt + afield
       ! Second and third steps are AB2 like
       afield = afield2%rhs() * (5. / 12.) * self%dt &
            & + afield%rhs() * (- 17. / 60.) * self%dt + afield2
       afield = afield%rhs() * (3. / 4.) * self%dt &
            & + afield2%rhs() * (- 5. / 12.) * self%dt + afield
    end do
  end subroutine integrate_RK3

  pure function euler_timestep(afield, dt) result(res)
    class(field_type), intent(in) :: afield
    real, intent(in) :: dt
    class(field_type), allocatable :: res
      res = afield%rhs() *  dt + afield
  end function euler_timestep

  subroutine AB2_timestep(f1, f2, dt)
    class(field_type), allocatable, intent(inout) :: f1, f2
    real, intent(in) :: dt

    f2 = f2%rhs() * 1.5 * dt - f1%rhs() * 0.5 * dt + f2
    f1 = f2
  end subroutine AB2_timestep

end module time_integrator
