program test_differentiate
  use iso_fortran_env, only: stderr => error_unit
  use differentiate, only: diff, diff2
  implicit none

  real :: f(20), df(20), expected(20)
  real, parameter :: tol = 0.1
  integer :: i, n
  logical :: allpass
  real :: dx

  n = size(f)
  dx = 2. * acos(-1.) / (n - 1)
  f = [(sin((i-1)*dx), i=1,n)]

  allpass = .true.

  ! First derivative
  expected = [(cos((i-1)*dx), i=1,n)]
  df = diff(f, dx)
  if (.not. all(abs(df - expected) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'First derivatives are computed correctly... failed'
  else
     write(stderr, '(a)') 'First derivatives are computed correctly... passed'
  end if

  ! Second derivative
  expected = [(- sin((i-1)*dx), i=1,n)]
  df = diff2(f, dx)
  if (.not. all(abs(df - expected) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'Second derivatives are computed correctly... failed'
  else
     write(stderr, '(a)') 'Second derivatives are computed correctly... passed'
  end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
  end if
end program test_differentiate
