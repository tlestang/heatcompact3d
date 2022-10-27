module differentiate
  use thomas_module, only: thomas
  implicit none

  integer, parameter :: dirichlet_nodes(8) = [ &
       & 1, 2, 3, 4, &
       & 1, 2, 3, 4 &
       & ]

  type :: differentiator_type
     private
     real :: coeffs(4, 2)
     real :: coeffs2(4, 2)
     integer :: nodes(4, 2)
   contains
     procedure, public :: diff, diff2
  end type differentiator_type

contains

  pure function dirichlet_differentiator()
    use coefficients_module, only: dirichlet_coeffs, &
         & dirichlet_2_coeffs
    type(differentiator_type) :: dirichlet_differentiator

    dirichlet_differentiator = differentiator_type( &
         & coeffs = reshape(dirichlet_coeffs, [4, 2]), &
         & coeffs2 = reshape(dirichlet_2_coeffs, [4, 2]), &
         & nodes = reshape(dirichlet_nodes, [4, 2]) &
         & )
  end function dirichlet_differentiator

  function neumann_odd_differentiator()
    use coefficients_module, only: neumann_odd_coeffs
    type(differentiator_type) :: neumann_odd_differentiator

    neumann_odd_differentiator = differentiator_type( &
         & coeffs = reshape(neumann_odd_coeffs, [4, 2]), &
         & coeffs2 = reshape(neumann_odd_coeffs, [4, 2]), &
         & nodes = reshape(dirichlet_nodes, [4, 2]) &
         & )
  end function neumann_odd_differentiator

  pure function diff(self, f, dx) result(df)
    use coefficients_module, only: sixth_order_compact_coeffs, alpha_coeffs
    class(differentiator_type), intent(in) :: self
    real, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: df(:), rhs(:), diag(:), as(:), cs(:)
    integer :: west_nodes(4, 2)
    integer :: nx, i
    real :: bulk_coeffs(4), east_coeffs(4, 2), west_coeffs(4, 2)

    allocate(rhs, source=f)

    east_coeffs = self%coeffs / dx
    rhs(1) = dot_product(f(self%nodes(:, 1)), east_coeffs(:, 1))
    rhs(2) = dot_product(f(self%nodes(:, 2)), east_coeffs(:, 2))

    nx = size(f)
    bulk_coeffs = sixth_order_compact_coeffs / dx
    do i=3,nx-2
       rhs(i) = dot_product(f([i-2, i-1, i+1, i+2]), bulk_coeffs)
    end do

    west_nodes = nx - self%nodes + 1
    west_coeffs = -1. * self%coeffs / dx
    rhs(nx - 1) = dot_product(f(west_nodes(:, 2)), west_coeffs(:, 2))
    rhs(nx) = dot_product(f(west_nodes(:, 1)), west_coeffs(:, 1))

    diag = [(1., i=1, nx)]
    cs = [2., 1. / 4., (alpha_coeffs, i=3, nx-2), 1. / 4., 0.]
    as = [0., 1. / 4., (alpha_coeffs, i=3, nx-2), 1. / 4., 2.]
    df = thomas(as, diag, cs, rhs)
  end function diff

  pure function diff2(self, f, dx) result(ddf)
    use coefficients_module, only: sixth_order_compact_2_coeffs, alpha_2_coeffs
    class(differentiator_type), intent(in) :: self
    real, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: ddf(:), rhs(:), diag(:), as(:), cs(:)
    integer :: west_nodes(4, 2)
    integer :: nx, i
    real :: bulk_coeffs(5), east_coeffs(4, 2), west_coeffs(4, 2)

    allocate(rhs, source=f)

    east_coeffs = self%coeffs2 / dx / dx
    rhs(1) = dot_product(f(self%nodes(:, 1)), east_coeffs(:, 1))
    rhs(2) = dot_product(f(self%nodes(:, 2)), east_coeffs(:, 2))

    nx = size(f)
    bulk_coeffs = sixth_order_compact_2_coeffs / dx / dx
    do i=3,nx-2
       rhs(i) = dot_product(f([i-2, i-1, i, i+1, i+2]), bulk_coeffs)
    end do

    west_nodes = nx - self%nodes + 1
    west_coeffs = + self%coeffs2 / dx / dx
    rhs(nx - 1) = dot_product(f(west_nodes(:, 2)), west_coeffs(:, 2))
    rhs(nx) = dot_product(f(west_nodes(:, 1)), west_coeffs(:, 1))

    diag = [(1., i=1, nx)]
    cs = [11., 1. / 10., (alpha_2_coeffs, i=3, nx-2), 1. / 10., 0.]
    as = [0., 1. / 10., (alpha_2_coeffs, i=3, nx-2), 1. / 10., 11.]
    ddf = thomas(as, diag, cs, rhs)
  end function diff2

end module differentiate
