program test_thomas
  use iso_fortran_env, only: stderr => error_unit
  use thomas_module, only: thomas
  implicit none

  logical :: allpass
  real, parameter :: tol = 0.0001
  real, parameter :: as(6) = [0., 1. / 4., 1. / 3., 1. / 3., 1. / 4., 2.]
  real, parameter :: cs(6) = [2., 1. / 4., 1. / 3., 1. / 3., 1. / 4., 0.]
  real, parameter :: bs(6) = [1., 1., 1., 1., 1., 1.]
  real, parameter :: ds(6) = [0., 1., 3., 3., 2., 1.]
  ! Solution computed with numpy.linalg.solve()
  real, parameter :: expected(6) = [ &
       & -1.71428, &
       & 0.85714, &
       & 2.28571, &
       & 1.28571, &
       & 2.85714, &
       & -4.71428 &
       & ]
  real :: x(6)

  allpass = .true.

  x = thomas(as, bs, cs, ds)
  if (.not. all(abs(x - expected) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Tridiagonal system is solved correctly... failed'
  else
     write(stderr, '(a)') 'Tridiagonal system is solved correctly... passed'
  end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
  end if

end program test_thomas
