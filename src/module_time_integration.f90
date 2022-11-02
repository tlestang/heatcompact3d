module time_integration
  use field_module, only: Field
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

  abstract interface
     subroutine integrate_proc(self, afield)
       import integrator_type, Field
       class(integrator_type), intent(in) :: self
       type(Field), intent(inout) :: afield
     end subroutine integrate_proc
  end interface

contains

  subroutine integrate_euler(self, afield)
    class(euler_integrator_type), intent(in) :: self
    type(Field), intent(inout) :: afield
    integer :: nt
    integer :: i

    nt = floor((self%endtime - self%starttime) / self%dt)
    do i = 1, nt
       afield = euler_timestep(afield, self%dt)
    end do
  end subroutine integrate_euler

  pure function euler_timestep(afield, dt) result(res)
    type(Field), intent(in) :: afield
    real, intent(in) :: dt
    type(Field) :: res
    res = afield%rhs() *  dt + afield
  end function euler_timestep
end module time_integration
