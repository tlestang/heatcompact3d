module thomas_module
  implicit none

contains

  pure function thomas(alpha, d) result(x)
    real, intent(in) :: alpha
    real, intent(in) :: d(:)
    real, allocatable :: b(:), x(:), dd(:)
    integer :: i, n
    real :: w

    allocate(dd, source=d)
    allocate(b, source=d)
    allocate(c, source=d)
    allocate(x, source=d)
    n = size(d)
    b = 1.
    do i = 2, n
       w = alpha / b(i-1)
       b(i) = 1. - w * alpha
       dd(i) = dd(i) - w * dd(i-1)
    end do

    x(n) = dd(n) / b(n)
    do i = n-1, 1, -1
       x(i) = (dd(i) - alpha * x(i+1)) / b(i)
    end do
  end function thomas
end module thomas_module
