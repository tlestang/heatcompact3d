program main
  use thomas_module, only: thomas
  implicit none

  real, parameter :: alpha = 2.
  real, parameter :: d(3) = [1, 2, 3]
  real :: x(3)

  x = thomas(alpha, d)

  ! Expect x = - 5. / 7.
  !        y =   6. / 7.
  !        z =   9. / 7.
  write(*,*) x
end program main
