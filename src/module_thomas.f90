module thomas_module
  implicit none

contains

  pure function thomas(as, bs, cs, ds) result(x)
    real, intent(in) :: as(:), bs(:), cs(:), ds(:)
    real, allocatable :: x(:), b(:), d(:)
    integer :: i, n
    real :: w

    allocate(d, source=ds)
    allocate(b, source=bs)
    allocate(x, source=ds)
    n = size(ds)
    do i = 2, n
       w = as(i) / b(i-1)
       b(i) = b(i) - w * cs(i-1)
       d(i) = d(i) - w * d(i-1)
    end do

    x(n) = d(n) / b(n)
    do i = n-1, 1, -1
       x(i) = (d(i) - cs(i) * x(i+1)) / b(i)
    end do
  end function thomas
end module thomas_module
