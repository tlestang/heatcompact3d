program test_stencil
  use iso_fortran_env, only: stderr => error_unit
  use stencil, only: stencil_type
  implicit none

  type(stencil_type) :: st, st_res, st_expected
  real, allocatable :: f(:)
  real :: coeffs(4)
  logical :: allpass
  real :: res_scal, expected_scal
  real, allocatable :: res(:), expected(:)

  real, parameter :: tol = 0.001

  coeffs = [2., 1., 2., 1.]
  st = stencil_type( &
       & nodes = [-1, 0, 1, 2], &
       & coeffs = coeffs, upper = 0., lower = 0. &
       & )
  f = [1., 2., 3., 4., 5., 6., 7., 8.]

  allpass = .true.

  res_scal = st%apply(f, ref=2)
  expected_scal = dot_product(f(1:4), coeffs)
  if (.not. (abs(res_scal - expected_scal) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Stencil is applied correctly... failed'
  else
     write(stderr, '(a)') 'Stencil is applied correctly... passed'
  end if

  res = st%apply_along(f)
  expected = [23., 14., 20., 26., 32., 38., 37., 29.]
  if (.not. all(abs(res - expected) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Stencil is applied along correctly... failed'
  else
     write(stderr, '(a)') 'Stencil is applied along correctly... passed'
  end if

  st_res = st%flip()
  st_expected = stencil_type( &
       & nodes = [+1, 0, -1, -2], &
       & coeffs = coeffs, upper = 0., lower = 0.&
       & )
  if (.not. st_expected%is_equal(st_res, tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Stencil is flipped correctly... failed'
  else
     write(stderr, '(a)') 'Stencil is flipped correctly... passed'
  end if

  st_res = st * 1.3
  st_expected = stencil_type( &
       & nodes = [-1, 0, 1, 2], &
       & coeffs = 1.3 * coeffs, upper = 0., lower = 0. &
       & )
  if (.not. st_expected%is_equal(st_res, tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Stencil is multiplied by scalar correctly... failed'
  else
     write(stderr, '(a)') 'Stencil is multiplied by scalar correctly... passed'
  end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
     error stop
  end if

end program test_stencil
