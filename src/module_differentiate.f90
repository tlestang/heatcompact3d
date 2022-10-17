module differentiate
  implicit none

  real, parameter :: afix = (7. / 9.)
  real, parameter :: bfix = (1. / 36.)

  real, parameter :: asix = 12. / 11.
  real, parameter :: bsix = 3. / 11.
  real, parameter :: csix = 0.

  real, parameter :: weights(4) = [ &
       - bfix, &
       & - afix, &
       & + afix, &
       & + bfix &
       & ]

  real, parameter :: weights2(7) = [ &
       & csix, &
       & bsix, &
       & asix, &
       & -2. * (asix + bsix + csix), &
       & asix, &
       & bsix, &
       & csix &
       & ]

contains

  pure function diff(f, dx) result(df)
    real, intent(in) :: f(:)
    real, intent(in) :: dx

    real, allocatable :: df(:)
    real, allocatable :: f_haloed(:)
    integer :: nx, i
    real :: w(4)

    nx = size(f)

    allocate(f_haloed(-1:nx+2))
    f_haloed(-1) = f(nx-2)
    f_haloed(0) = f(nx-1)
    f_haloed(nx + 1) = f(2)
    f_haloed(nx + 2) = f(3)
    f_haloed(1:nx) = f

    w = weights / dx
    allocate(df, source=f)
    do i=1,nx
       df(i) = sum(f_haloed([i-2, i-1, i+1, i+2]) * w)
    end do
  end function diff


  pure function diff2(f, dx) result(ddf)
    real, intent(in) :: f(:)
    real, intent(in) :: dx

    real, allocatable :: ddf(:)
    real, allocatable :: f_haloed(:)
    integer :: nx, i
    real :: w(7)

    nx = size(f)

    allocate(f_haloed(-2:nx+3))
    f_haloed(-2) = f(nx-3)
    f_haloed(-1) = f(nx-2)
    f_haloed(0) = f(nx-1)
    f_haloed(nx + 1) = f(2)
    f_haloed(nx + 2) = f(3)
    f_haloed(nx + 3) = f(4)

    w = weights2 / dx / dx
    allocate(ddf, source=f)
    do i=1,nx
       ddf(i) = sum(f_haloed([i-3, i-2, i-1, i, i+1, i+2, i+3]) * w)
    end do
  end function diff2

end module differentiate
