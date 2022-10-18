program main
  use thomas_module, only: thomas
  implicit none

  real, parameter :: alpha = 2.
  real, parameter :: ds(3) = [1, 2, 3]
  real, parameter :: bs(3) = [1, 1, 1]
  real :: x(3)

  x = thomas(bs, ds, alpha)

  ! Expect x = - 5. / 7.
  !        y =   6. / 7.
  !        z =   9. / 7.
  write(*,*) x
end program main
