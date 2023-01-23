module differentiate
  !! This module implements differentiation over a one dimensional
  !! pencil.  It exports concrete derived types
  !!
  !! - `differentiator_type` :: periodic differentiation
  !! - `nonperiodic_differentiator_type` :: non periodic differentiation
  !!
  !! Both types provide access to a single type-procedure `diff` that
  !! performs the differentiation and return the derivative evaluated
  !! over the input pencil.
  !!
  !! Instances of differentiator derived types are created on the
  !! client side via constructor procedures `<method>_<order>`
  !! (e.g. `sixth_order_compact_1`). These constructor procedures are
  !! implemented as a generic interface to concrete constructors for
  !! the periodic case (no boundary conditions are specified) or the
  !! noneriodic case (boundary conditions are specified). See module
  !! `boundary_schemes`.
  !!
  !! ```
  !! use boundary_schemes, only:: get_dirichlet_4th
  !! real :: f(:), df(:), dx
  !! ! ...
  !! ! ...
  !! differentiatior = sixth_order_compact_1( &
  !!      & east = get_dirichlet_4th(), &
  !!      & west = get_dirichlet_4th(), &
  !!      & )
  !! df = differentiatior%diff(f, dx)
  !! ```

  use stencil, only: stencil_type
  use boundary_schemes, only: boundary_type, sixth_order_compact_stencil, &
       & sixth_order_compact_second_stencil
  use thomas_module, only: thomas
  implicit none

  type :: differentiator_type
     !! Implements differentiation over a periodic stencil.
     !! `differentiator_type` provides access a a unique type-bound
     !! procedure `diff`
     private
     type(stencil_type) :: bulk_stencil !! Differentiation stencil
   contains
     procedure, public :: diff => diff_periodic
  end type differentiator_type

  type, extends(differentiator_type) :: nonperiodic_differentiator_type
     type(stencil_type) :: east_stencils(2)
     type(stencil_type) :: west_stencils(2)
   contains
     procedure, public :: diff => diff_nonperiodic
  end type nonperiodic_differentiator_type

  interface sixth_order_compact_1
     module procedure :: sixth_order_compact_1_periodic
     module procedure :: sixth_order_compact_1_nonperiodic
  end interface sixth_order_compact_1

  interface sixth_order_compact_2
     module procedure :: sixth_order_compact_2_periodic
     module procedure :: sixth_order_compact_2_nonperiodic
  end interface sixth_order_compact_2

contains

  pure function sixth_order_compact_1_nonperiodic(east, west)
    type(boundary_type), intent(in) :: east, west
    type(nonperiodic_differentiator_type) :: sixth_order_compact_1_nonperiodic
    sixth_order_compact_1_nonperiodic = nonperiodic_differentiator_type( &
         & east_stencils = east%first_order_east, &
         & west_stencils = west%first_order_west, &
         & bulk_stencil = sixth_order_compact_stencil() &
         & )
  end function sixth_order_compact_1_nonperiodic

  pure function sixth_order_compact_1_periodic()
    type(differentiator_type) :: sixth_order_compact_1_periodic
    sixth_order_compact_1_periodic = differentiator_type( &
         & bulk_stencil = sixth_order_compact_stencil() &
         & )
  end function sixth_order_compact_1_periodic

  pure function sixth_order_compact_2_nonperiodic(east, west)
    type(boundary_type), intent(in) :: east, west
    type(nonperiodic_differentiator_type) :: sixth_order_compact_2_nonperiodic
    sixth_order_compact_2_nonperiodic = nonperiodic_differentiator_type( &
         & east_stencils = east%second_order_east, &
         & west_stencils = west%second_order_west, &
         & bulk_stencil = sixth_order_compact_second_stencil() &
         & )
  end function sixth_order_compact_2_nonperiodic

  pure function sixth_order_compact_2_periodic()
    type(differentiator_type) :: sixth_order_compact_2_periodic
    sixth_order_compact_2_periodic = differentiator_type( &
         & bulk_stencil = sixth_order_compact_second_stencil() &
         & )
  end function sixth_order_compact_2_periodic

  pure function diff_nonperiodic(self, f, dx) result(df)
    !! Apply a differentiation stencil along a one dimensional pencil,
    !! then apply boundary conditions on both ends. Boundary
    !! conditions are applied as arrays of type `stencil_type`

    class(nonperiodic_differentiator_type), intent(in) :: self
    real, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: df(:), rhs(:), diag(:), lower_diag(:), upper_diag(:)
    integer :: neast, nwest
    integer :: ref, i ! do loop counters
    type(stencil_type) :: sten
    real :: alpha_coeffs

    neast = size(self%east_stencils)
    nwest = size(self%west_stencils)

    ! Compute rhs assuming periodicity
    rhs = self%bulk_stencil%apply_along(f)

    ! Fix rhs values near boundaries
    do ref = 1, neast
       sten = self%east_stencils(ref)
       rhs(ref) = sten%apply(f, ref)
    end do

    do i = 1, nwest
       sten = self%west_stencils(i)
       ref = size(f) - i + 1
       rhs(ref) = sten%apply(f, ref)
    end do

    rhs = rhs / dx
    ! Solve tridiagonal system of equations
    diag = [(1., i=1, size(f))]
    ! Both upper and lower diagonals are declared of size n = size(f)
    ! instead of n-1, out of convenience. Components upper_diag(n) and
    ! lower_diag(1) will not be accessed by the thomas solver as they
    ! do not appear in the tridiagonal system.
    upper_diag = [ &
         & self%east_stencils%get_upper(), &
         & (self%bulk_stencil%get_upper(), i=neast + 1, size(f) - nwest), &
         & reverse(self%west_stencils%get_upper()) &
         ]
    lower_diag = [ &
         & self%east_stencils%get_upper(), &
         & (self%bulk_stencil%get_lower(), i=neast + 1, size(f) - nwest), &
         & reverse(self%west_stencils%get_upper()) &
         ]
    df = thomas(lower_diag, diag, upper_diag, rhs)
  end function diff_nonperiodic

  pure function diff_periodic(self, f, dx) result(df)
    !! Apply a differentiation stencil along a one dimensional pencil,
    !! assuming periodic boundaries. For instance, with a four point
    !! stencil \(s = (-2, -1, 1, 2)\) and weights \({g_{i}}_{1 \leq i
    !! \leq 4}\)
    !!
    !! \[ g(1) = g_1f_{n-2} + g_2f_{n-1} + g_3f_{2} + g_4f_{3} \]
    !! \[ g(2) = g_1f_{n-1} + g_2f_{1} + g_3f_{3} + g_4f_{4} \]
    !!
    !! \[ g(i) = g_1f_{i-2} + g_2f_{i-1} + g_3f_{i+1} + g_4f_{i+2}, i = 3, n-2\]
    !!
    !! \[ g(n-1) = g_1f_{n-3} + g_2f_{n-2} + g_3f_{n} + g_4f_{2} \]
    !! \[ g(n) = g_1f_{n-2} + g_2f_{n-1} + g_3f_{2} + g_4f_{3} \]
    class(differentiator_type), intent(in) :: self
    real, intent(in) :: f(:)
    !! Function to be derive, evaluated on pencil
    real, allocatable :: df(:)
    !! Derivative, evaluated on pencil
    real, intent(in) :: dx
    !! Step size

    real, allocatable :: rhs(:), diag(:), lower_diag(:), upper_diag(:)
    real, allocatable :: u(:), v(:), q(:), y(:)
    real :: gamma, up, low
    integer :: nx, i

    rhs = self%bulk_stencil%apply_along(f)
    rhs = rhs / dx
    ! Solve quasi-tridiagonal system of equations using the
    ! Shermann-Morrison formula
    up = self%bulk_stencil%get_upper()
    low = self%bulk_stencil%get_lower()
    gamma = -1.
    nx = size(f)
    diag = [ &
         & 1. - gamma, &
         & (1., i=2, nx - 1), &
         & 1. - (up * low) / gamma &
         & ]
    ! Both upper and lower diagonals are declared of size n = size(f)
    ! instead of n-1, out of convenience. Components upper_diag(n) and
    ! lower_diag(1) will not be accessed by the thomas solver as they
    ! do not appear in the tridiagonal system.
    nx = size(f)

    lower_diag = [(low, i = 1, nx)]
    upper_diag = [(up, i = 1, nx)]
    u = [gamma, (0., i=2, nx-1), up]
    v = [1., (0., i=2, nx-1), low / gamma]
    q = thomas(lower_diag, diag, upper_diag, u)
    y = thomas(lower_diag, diag, upper_diag, rhs)

    df = y - ((y(1) - low * y(nx)) &
         & / (1. + q(1) - low * q(nx))) * q
  end function diff_periodic

  pure function reverse(x)
    real, intent(in) :: x(:)
    real :: reverse(size(x))
    reverse = x(size(x):1:-1)
  end function reverse

end module differentiate
