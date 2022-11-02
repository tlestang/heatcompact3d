program test_field
  use iso_fortran_env, only: stderr => error_unit
  use field_module, only: Field
  use time_integration, only: euler_integrator_type
  implicit none

  real :: u0(16, 16, 16)
  type(Field) :: temp_field, expected
  real, parameter :: tol = 0.1
  integer :: i, j, k, nx, ny, nz
  logical :: allpass
  real :: dx
  type(euler_integrator_type) :: euler

  nx = size(u0, 1)
  ny = size(u0, 2)
  nz = size(u0, 3)
  dx = 2. * acos(-1.) / (nx - 1)

  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           u0(i, j, k) = &
                & cos((i-1)*dx) * cos((j-1)*dx) * &
                & cos((k-1)*dx)
        end do
     end do
  end do

  allpass = .true.
  dt = 0.1

  temp_field = Field(u0, dx)
  euler = euler_integrator_type(starttime=0., endtime=0.3, dt=dt)
  call euler%integrate(temp_field)
  expected = Field( &
       & (1. - 3. * dt) ** 3 * u0, dx)
  if (.not. expected%is_equal(temp_field, tol)) then
     write(stderr, '(a)') 'Foward integration (Euler) is computed correctly... failed.'
     allpass = .false.
  else
     write(stderr, '(a)') 'Foward integration (Euler) is computed correctly... passed.'
  end if
  end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
  end if

end program test_field
