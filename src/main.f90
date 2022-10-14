program main
  use differentiate, only: diff
  implicit none

  real :: f(100), df(100)
  integer, parameter :: nx = 100
  integer :: i
  real :: dx

  dx = 2. * 3.14 / (nx - 1)
  f = [(sin((i-1)*dx), i=1,nx)]
  df = diff(f, dx)
  write(*,*) df

end program main
