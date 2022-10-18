program main
  use differentiate, only: diff2, diff
  implicit none

  real :: f(100), ddf(100), expected_ddf(100)
  integer, parameter :: nx = 100
  integer :: i
  real :: dx

  dx = 2. * acos(-1.) / (nx - 1)
  f = [(sin((i-1)*dx), i=1,nx)] 
  expected_ddf = [(- sin((i-1)*dx), i=1,nx)]
  ddf = diff2(f, dx)
  write(*,*) ddf
  write(*,*) expected_ddf

end program main
