program test_differentiate
  use iso_fortran_env, only: stderr => error_unit
  use differentiate, only: differentiator_type, &
       & sixth_order_compact
  use boundary_schemes, only: boundary_type, get_dirichlet_boundary
  implicit none

  real, allocatable :: f(:), df(:), expected(:)
  real, parameter :: tol = 0.1
  type(differentiator_type) :: differentiator
  type(boundary_type) :: dirichlet
  integer :: i, n
  logical :: allpass
  real :: dx

  n = 20
  dx = 2. * acos(-1.) / (n - 1)
  f = [(sin((i-1)*dx), i=1,n)]

  allpass = .true.

  ! First derivative with dirichlet both ends
  expected = [(cos((i-1)*dx), i=1,n)]
  dirichlet = get_dirichlet_boundary()
  differentiator = sixth_order_compact( &
       & east = dirichlet, &
       & west = dirichlet &
       & )
  df = differentiator%diff(f, dx)
  if (.not. all(abs(df - expected) < tol)) then
     allpass = .false.
     write(stderr, '(a)') 'First derivatives are computed correctly... failed'
  else
     write(stderr, '(a)') 'First derivatives are computed correctly... passed'
  end if

  ! ! Second derivative
  ! expected = [(- sin((i-1)*dx), i=1, n)]
  ! df = differentiator%diff2(f, dx)
  ! if (.not. all(abs(df - expected) < tol)) then
  !    allpass = .false.
  !    write(stderr, '(a)') 'Second derivatives are computed correctly... failed'
  ! else
  !    write(stderr, '(a)') 'Second derivatives are computed correctly... passed'
  ! end if

  ! dx = acos(-1.) / (n - 1)
  ! f = [(cos((i-1)*dx), i=1,n)]
  ! ! First derivative with neumann odd both ends
  ! expected = [(-sin((i-1)*dx), i=1,n)]
  ! differentiator = neumann_odd_differentiator()
  ! df = differentiator%diff(f, dx)
  ! if (.not. all(abs(df - expected) < tol)) then
  !    allpass = .false.
  !    write(stderr, '(a)') 'First derivatives are computed correctly... failed'
  ! else
  !    write(stderr, '(a)') 'First derivatives are computed correctly... passed'
  ! end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
  end if
end program test_differentiate
