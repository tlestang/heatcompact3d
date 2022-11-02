program test_field
  use iso_fortran_env, only: stderr => error_unit
  use field_module, only: Field
  use time_integration, only: euler_integrator_type, AB2_integrator_type, &
       & RK3_integrator_type
  implicit none

  real :: u0(16, 16, 16)
  type(Field) :: temp_field, expected
  real, parameter :: tol = 0.1
  integer :: i, j, k, nx, ny, nz
  logical :: allpass
  real :: dx, dt, dt2, dt3
  type(euler_integrator_type) :: euler
  type(AB2_integrator_type) :: AB2
  type(RK3_integrator_type) :: RK3

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
  dt2 = dt * dt
  dt3 = dt * dt * dt

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

  temp_field = Field(u0, dx)
  AB2 = AB2_integrator_type(starttime=0., endtime=0.3, dt=dt)
  call AB2%integrate(temp_field)
  expected = Field( &
       & 0.25 * (-243.*dt3 + 144.*dt2 - 36.*dt + 4) * u0, dx)
  if (.not. expected%is_equal(temp_field, tol)) then
     write(stderr, '(a)') 'Foward integration (Adams-Bashforth 2) is computed correctly... failed.'
     allpass = .false.
  else
     write(stderr, '(a)') 'Foward integration (Adams-Bashforth 2) is computed correctly... passed.'
  end if

  temp_field = Field(u0, dx)
  RK3 = RK3_integrator_type(starttime=0., endtime=0.2, dt=dt)
  call RK3%integrate(temp_field)
  expected = Field( &
       &  (9.*dt3 - 9.*dt2 + 6.*dt - 2)**2 * u0 / 4., dx)
  if (.not. expected%is_equal(temp_field, tol)) then
     write(stderr, '(a)') 'Foward integration (Runge-Kutta 3) is computed correctly... failed.'
     allpass = .false.
  else
     write(stderr, '(a)') 'Foward integration (Runge-Kutta 3) is computed correctly... passed.'
  end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
  end if

end program test_field
