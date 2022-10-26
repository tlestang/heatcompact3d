module differentiate
  use thomas_module, only: thomas
  implicit none

  real, parameter :: dirichlet_coeffs(8) = [ &
       & -5. / 2., 2., 0.5, 0., &
       & -3. / 4., 0., 3. / 4., 0. &
       & ]
  integer, parameter :: dirichlet_stencils(8) = [ &
       & 1, 2, 3, 4, &
       & 1, 2, 3, 4 &
       & ]
  real, parameter :: afix = (7. / 9.)
  real, parameter :: bfix = (1. / 36.)
  real, parameter :: sixth_order_compact(4) = [ &
       - bfix, &
       & - afix, &
       & + afix, &
       & + bfix &
       & ]
  real, parameter :: alpha = 1. / 3.

  real, parameter :: asix = 12. / 11.
  real, parameter :: bsix = 3. / 44.
  real, parameter :: csix = 0.
  real, parameter :: alpha2 = 2. / 11.

  real, parameter :: weights2(7) = [ &
       & csix, &
       & bsix, &
       & asix, &
       & -2. * (asix + bsix + csix), &
       & asix, &
       & bsix, &
       & csix &
       & ]

  real, parameter :: gamma = -1.

  type :: differentiator_type
     private
     real :: east_coeffs(4, 2)
     real :: west_coeffs(4, 2)
     integer :: east_stencils(4, 2)
   contains
     procedure, public :: diff
  end type differentiator_type

contains

  function dirichlet_differentiator()
    type(differentiator_type) :: dirichlet_differentiator

    dirichlet_differentiator = differentiator_type( &
         & east_coeffs = reshape(dirichlet_coeffs, [4, 2]), &
         & west_coeffs = reshape(dirichlet_coeffs, [4, 2]), &
         & east_stencils = reshape(dirichlet_stencils, [4, 2]) &
         & )
  end function dirichlet_differentiator

  pure function diff(self, f, dx) result(df)
    class(differentiator_type), intent(in) :: self
    real, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: df(:), rhs(:), diag(:), as(:), cs(:)
    integer :: west_stencils(4, 2)
    integer :: nx, i
    real :: w(4)

    allocate(rhs, source=f)

    rhs(1) = dot_product(f(self%east_stencils(:, 1)), self%east_coeffs(:, 1) / dx)
    rhs(2) = dot_product(f(self%east_stencils(:, 2)), self%east_coeffs(:, 2) / dx)

    nx = size(f)
    w = sixth_order_compact / dx
    do i=3,nx-2
       rhs(i) = dot_product(f([i-2, i-1, i+1, i+2]), w)
    end do

    west_stencils = nx - self%east_stencils + 1
    rhs(nx - 1) = dot_product(f(west_stencils(:, 2)), - self%east_coeffs(:, 2) / dx)
    rhs(nx) = dot_product(f(west_stencils(:, 1)), - self%east_coeffs(:, 1) / dx)

    diag = [(1., i=1, nx)]
    cs = [2., 1. / 4., (alpha, i=3, nx-2), 1. / 4., 0.]
    as = [0., 1. / 4., (alpha, i=3, nx-2), 1. / 4., 2.]
    df = thomas(as, diag, cs, rhs)
  end function diff


  pure function diff2(f, dx) result(ddf)
    real, intent(in) :: f(:)
    real, intent(in) :: dx

    real, allocatable :: ddf(:)
    real, allocatable :: rhs(:)
    real, allocatable :: f_haloed(:)
    real, allocatable :: diag(:), u(:), v(:), q(:), y(:)
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
    f_haloed(1:nx) = f

    w = weights2 / (dx * dx)
    allocate(rhs, source=f)
    do i=1,nx
       rhs(i) = sum(f_haloed([i-3, i-2, i-1, i, i+1, i+2, i+3]) * w)
    end do
    diag = [1. - gamma, (1., i=2, nx - 1), 1. + alpha2 * alpha2]
    u = [gamma, (0., i=2, nx-1), alpha2]
    v = [1., (0., i=2, nx-1), - alpha2]
    q = thomas(diag, u, alpha2)
    y = thomas(diag, rhs, alpha2)

    ddf = y - ((y(1) - alpha2 * y(nx)) / (1. + q(1) - alpha2 * q(nx))) * q
  end function diff2

end module differentiate
