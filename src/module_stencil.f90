module stencil
  implicit none
  type :: stencil_type
     integer :: nodes(4)
     real :: coeffs(4)
   contains
     procedure, private :: stencil_mul_real
     procedure, public :: apply => apply_stencil
     procedure, public :: apply_along => apply_stencil_along
     procedure, public :: flip, is_equal
     generic :: operator(*) => stencil_mul_real
  end type stencil_type

  type(stencil_type), parameter :: sixth_order_compact_stencil = &
       stencil_type( &
       & nodes = [-2, -1, 1, 2], &
       & coeffs = [- 1. / 36., - 7. / 9., + 7. / 9., + 1. / 36.] &
       & )

contains

  pure logical function is_equal(self, st, tol)
    class(stencil_type), intent(in) :: self
    type(stencil_type), intent(in) :: st
    real, intent(in) :: tol
    logical nodes_equal, coeffs_equal
    nodes_equal = all(self%nodes == st%nodes)
    coeffs_equal = all(abs(self%coeffs - st%coeffs) < tol)

    is_equal = nodes_equal .and. coeffs_equal
  end function is_equal

  pure elemental type(stencil_type) function stencil_mul_real(self, a)
    class(stencil_type), intent(in) :: self
    real, intent(in) :: a
    stencil_mul_real = stencil_type( &
         & nodes = self%nodes, &
         & coeffs = a * self%coeffs &
         & )
  end function stencil_mul_real

  pure elemental type(stencil_type) function flip(self)
    class(stencil_type), intent(in) :: self
    flip = stencil_type( &
         & nodes = - self%nodes, &
         & coeffs = self%coeffs &
         & )
  end function flip

  pure real function apply_stencil(self, f, ref)
    class(stencil_type), intent(in) :: self
    real, allocatable, intent(in) :: f(:)
    integer, intent(in) :: ref
    real, allocatable :: eval(:)

    eval = f(self%nodes + ref)
    apply_stencil = dot_product(eval, self%coeffs)
  end function apply_stencil

  pure function apply_stencil_along(self, f)
    class(stencil_type), intent(in) :: self
    real, allocatable, intent(in) :: f(:)
    real, allocatable :: apply_stencil_along(:), f_padded(:)
    integer :: ref, i, lpad, rpad

    lpad = minval(self%nodes) + 1
    rpad = size(f) + maxval(self%nodes)
    allocate(f_padded(lpad:rpad))
    do i = lpad, 0
       f_padded(i) = f(size(f) + i)
    end do
    do i = size(f) + 1, rpad
       f_padded(i) = f(i - size(f))
    end do
    f_padded(1:size(f)) = f
    allocate(apply_stencil_along, source = f)
    do ref = 1, size(f)
       apply_stencil_along(ref) = self%apply(f_padded, ref)
    end do
  end function apply_stencil_along

end module stencil

module boundary_schemes
  use stencil, only: stencil_type
  type boundary_type
     type(stencil_type) :: first_order_east(2), second_order_east(2)
     type(stencil_type) :: first_order_west(2), second_order_west(2)
  end type boundary_type

contains

  type(boundary_type) function get_dirichlet_boundary()
    type(stencil_type) :: first_order_east(2), second_order_east(2)
    type(stencil_type) :: first_order_west(2), second_order_west(2)

    first_order_east(1) = stencil_type( &
         & nodes = [0, 1, 2, 3], &
         & coeffs = [-5. / 2., 2., 0.5, 0.] &
         & )
    first_order_east(2) = stencil_type( &
         & nodes = [-1, 0, 1, 2], &
         & coeffs = [-3. / 4., 0., 3. / 4., 0.] &
         & )
    second_order_east(1) = stencil_type( &
         & nodes = [0, 1, 2, 3], &
         & coeffs = [13., -27., 15., -1.] &
         & )
    second_order_east(2) = stencil_type( &
         & nodes = [-1, 0, 1, 2], &
         & coeffs = [6. / 5., -12. / 5., 6. / 5., 0.] &
         & )
    first_order_west = first_order_east%flip() * (-1.)
    second_order_west = second_order_east%flip()
    get_dirichlet_boundary = boundary_type( &
         & first_order_east, second_order_east, first_order_west, second_order_west)
  end function get_dirichlet_boundary
end module boundary_schemes
