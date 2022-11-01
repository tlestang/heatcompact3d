module differentiate
  use stencil, only: stencil_type
  use boundary_schemes, only: boundary_type, sixth_order_compact_stencil, &
       & sixth_order_compact_second_stencil
  use thomas_module, only: thomas
  implicit none

  type :: differentiator_type
     private
     type(stencil_type) :: bulk_stencil
   contains
     procedure, public :: diff => diff_periodic
  end type differentiator_type

  type, extends(differentiator_type) :: nonperiodic_differentiator_type
     type(stencil_type) :: east_stencils(2)
     type(stencil_type) :: west_stencils(2)
   contains
     procedure, public :: diff => diff_nonperiodic
  end type nonperiodic_differentiator_type

  interface sixth_order_compact
     module procedure :: sixth_order_compact_periodic
     module procedure :: sixth_order_compact_nonperiodic
  end interface sixth_order_compact

  interface sixth_order_compact_2nd
     module procedure :: sixth_order_compact_periodic_2nd
     module procedure :: sixth_order_compact_nonperiodic_2nd
  end interface sixth_order_compact_2nd

contains

  pure function sixth_order_compact_nonperiodic(east, west)
    type(boundary_type), intent(in) :: east, west
    type(nonperiodic_differentiator_type) :: sixth_order_compact_nonperiodic
    sixth_order_compact_nonperiodic = nonperiodic_differentiator_type( &
         & east_stencils = east%first_order_east, &
         & west_stencils = west%first_order_west, &
         & bulk_stencil = sixth_order_compact_stencil() &
         & )
  end function sixth_order_compact_nonperiodic

  pure function sixth_order_compact_periodic()
    type(differentiator_type) :: sixth_order_compact_periodic
    sixth_order_compact_periodic = differentiator_type( &
         & bulk_stencil = sixth_order_compact_stencil() &
         & )
  end function sixth_order_compact_periodic

  pure function sixth_order_compact_nonperiodic_2nd(east, west)
    type(boundary_type), intent(in) :: east, west
    type(nonperiodic_differentiator_type) :: sixth_order_compact_nonperiodic_2nd
    sixth_order_compact_nonperiodic_2nd = nonperiodic_differentiator_type( &
         & east_stencils = east%second_order_east, &
         & west_stencils = west%second_order_west, &
         & bulk_stencil = sixth_order_compact_second_stencil() &
         & )
  end function sixth_order_compact_nonperiodic_2nd

  pure function sixth_order_compact_periodic_2nd()
    type(differentiator_type) :: sixth_order_compact_periodic_2nd
    sixth_order_compact_periodic_2nd = differentiator_type( &
         & bulk_stencil = sixth_order_compact_second_stencil() &
         & )
  end function sixth_order_compact_periodic_2nd

  pure function diff_nonperiodic(self, f, dx) result(df)
    class(nonperiodic_differentiator_type), intent(in) :: self
    real, allocatable, intent(in) :: f(:)
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
    class(differentiator_type), intent(in) :: self
    real, allocatable, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: df(:), rhs(:), diag(:), lower_diag(:), upper_diag(:)
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
    real, allocatable :: reverse(:)
    reverse = x(size(x):1:-1)
  end function reverse

end module differentiate
