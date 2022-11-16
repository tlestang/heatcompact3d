module time_integration
  use field_module, only: field_type
  implicit none

  type, abstract :: integrator_type
     real :: starttime, endtime, dt
   contains
     procedure(integrate_proc), public, &
          & deferred :: integrate
  end type integrator_type

  type, extends(integrator_type) :: euler_integrator_type
     real :: alpha = 1.
   contains
     procedure :: integrate => integrate_euler
  end type euler_integrator_type

  type, extends(integrator_type) :: AB2_integrator_type
   contains
     procedure :: integrate => integrate_AB2
  end type AB2_integrator_type

  type, extends(integrator_type) :: RK3_integrator_type
   contains
     procedure :: integrate => integrate_RK3
  end type RK3_integrator_type

  abstract interface
     subroutine integrate_proc(self, afield)
       import integrator_type, field_type
       class(integrator_type), intent(in) :: self
       type(field_type), intent(inout) :: afield
     end subroutine integrate_proc
  end interface

contains

  subroutine integrate_euler(self, afield)
    class(euler_integrator_type), intent(in) :: self
    type(field_type), intent(inout) :: afield
    integer :: nt
    integer :: i

    nt = floor((self%endtime - self%starttime) / self%dt)
    do i = 1, nt
       afield = euler_timestep(afield, self%dt)
    end do
  end subroutine integrate_euler

  subroutine integrate_AB2(self, afield)
    class(AB2_integrator_type), intent(in) :: self
    type(field_type), intent(inout) :: afield
    type(field_type) :: fields(2)
    integer :: nt
    integer :: i

    nt = floor((self%endtime - self%starttime) / self%dt)
    fields(1) = afield
    fields(2) = euler_timestep(afield, self%dt)
    do i = 2, nt
       fields = AB2_timestep(fields, self%dt)
    end do
    afield = fields(2)
  end subroutine integrate_AB2

  subroutine integrate_RK3(self, afield)
    class(RK3_integrator_type), intent(in) :: self
    type(field_type), intent(inout) :: afield
    type(field_type) :: afield2
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
    type(field_type), intent(in) :: afield
    real, intent(in) :: dt
    type(field_type) :: res
    res = afield%rhs() *  dt + afield
  end function euler_timestep

  pure function AB2_timestep(fields, dt) result(res)
    type(field_type), intent(in) :: fields(2)
    real, intent(in) :: dt
    type(field_type) :: res(2)
    type(field_type) :: f1, f2

    f1 = fields(1)
    f2 = fields(2)
    res(2) = f2%rhs() * 1.5 * dt - f1%rhs() * 0.5 * dt + f2
    res(1) = f2
  end function AB2_timestep

end module time_integration
