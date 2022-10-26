module differentiate
  use thomas_module, only: thomas
  implicit none

  real, parameter :: dirichlet_coeffs(8) = [ &
       & -5. / 2., 2., 0.5, 0., &
       & -3. / 4., 0., 3. / 4., 0. &
       & ]
  real, parameter :: dirichlet_coeffs2(8) = [ &
       & 13., -27., 15., -1., &
       & 6. / 5., -12. / 5., 6. / 5., 0. &
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
  real, parameter :: neumann_odd_coeffs(8) = [ &
       & 0., 0., 0., 0., &
       & -afix, -bfix, afix, bfix &
       & ]
  real, parameter :: alpha = 1. / 3.

  real, parameter :: asix = 12. / 11.
  real, parameter :: bsix = 3. / 44.
  real, parameter :: sixth_order_compact2(5) = [ &
       & bsix, &
       & asix, &
       & -2. * (asix + bsix), &
       & asix, &
       & bsix &
       & ]
  real, parameter :: alpha2 = 2. / 11.
  real, parameter :: gamma = -1.

  type :: differentiator_type
     private
     real :: east_coeffs(4, 2)
     real :: east_coeffs2(4, 2)
     real :: west_coeffs(4, 2)
     integer :: east_stencils(4, 2)
   contains
     procedure, public :: diff, diff2
  end type differentiator_type

contains

  function dirichlet_differentiator()
    type(differentiator_type) :: dirichlet_differentiator

    dirichlet_differentiator = differentiator_type( &
         & east_coeffs = reshape(dirichlet_coeffs, [4, 2]), &
         & east_coeffs2 = reshape(dirichlet_coeffs2, [4, 2]), &
         & west_coeffs = reshape(dirichlet_coeffs, [4, 2]), &
         & east_stencils = reshape(dirichlet_stencils, [4, 2]) &
         & )
  end function dirichlet_differentiator

  function neumann_odd_differentiator()
    type(differentiator_type) :: neumann_odd_differentiator

    neumann_odd_differentiator = differentiator_type( &
         & east_coeffs = reshape(neumann_odd_coeffs, [4, 2]), &
         & east_coeffs2 = reshape(neumann_odd_coeffs, [4, 2]), &
         & west_coeffs = reshape(neumann_odd_coeffs, [4, 2]), &
         & east_stencils = reshape(dirichlet_stencils, [4, 2]) &
         & )
  end function neumann_odd_differentiator

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

  pure function diff2(self, f, dx) result(ddf)
    class(differentiator_type), intent(in) :: self
    real, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: ddf(:), rhs(:), diag(:), as(:), cs(:)
    integer :: west_stencils(4, 2)
    integer :: nx, i
    real :: w(5)

    allocate(rhs, source=f)

    rhs(1) = dot_product(f(self%east_stencils(:, 1)), self%east_coeffs2(:, 1) / dx / dx)
    rhs(2) = dot_product(f(self%east_stencils(:, 2)), self%east_coeffs2(:, 2) / dx / dx)

    nx = size(f)
    w = sixth_order_compact2 / dx / dx
    do i=3,nx-2
       rhs(i) = dot_product(f([i-2, i-1, i, i+1, i+2]), w)
    end do

    west_stencils = nx - self%east_stencils + 1
    rhs(nx - 1) = dot_product(f(west_stencils(:, 2)), self%east_coeffs2(:, 2) / dx / dx)
    rhs(nx) = dot_product(f(west_stencils(:, 1)), self%east_coeffs2(:, 1) / dx / dx)

    diag = [(1., i=1, nx)]
    cs = [11., 1. / 10., (alpha2, i=3, nx-2), 1. / 10., 0.]
    as = [0., 1. / 10., (alpha2, i=3, nx-2), 1. / 10., 11.]
    ddf = thomas(as, diag, cs, rhs)
  end function diff2

end module differentiate
