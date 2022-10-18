module thomas_module
  implicit none

contains

  pure function thomas(alpha, d) result(x)
    real, allocatable, intent(inout) :: d(:),
    real, allocatable :: b(;), c(:), x(:)
    do i = 2, n
       b(i) = 1. - alpha * alpha
       d(i) = d(i) - alpha * d(i-1)
    end do

    x(n) = d(n)
    do i = n-1, 1
       x(i) = d(i) - c(i) * x(i+1)
    end do
  end function thomas
end module thomas_module
