module differentiate
  use stencil, only: stencil_type, sixth_order_compact_stencil
  use boundary_schemes, only: boundary_type
  use thomas_module, only: thomas
  implicit none

  type :: differentiator_type
     private
     type(stencil_type) :: east_stencils(2)
     type(stencil_type) :: west_stencils(2)
     type(stencil_type) :: bulk_stencil
   contains
     procedure, public :: diff
  end type differentiator_type

  type, extends(differentiator_type) :: periodic_differentiator_type
   contains
     procedure, public :: diff => diff_periodic
  end type periodic_differentiator_type

contains

  pure function sixth_order_compact(east, west)
    type(boundary_type), intent(in) :: east, west
    type(differentiator_type) :: sixth_order_compact

    sixth_order_compact = differentiator_type( &
         & east_stencils = east%first_order_east, &
         & west_stencils = west%first_order_west, &
         & bulk_stencil = sixth_order_compact_stencil &
         & )
  end function sixth_order_compact

  pure function diff(self, f, dx) result(df)
    class(differentiator_type), intent(in) :: self
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
  end function diff

  pure function diff_periodic(self, f, dx) result(df)
    class(periodic_differentiator_type), intent(in) :: self
    real, allocatable, intent(in) :: f(:)
    real, intent(in) :: dx
    real, allocatable :: df(:), rhs(:), diag(:), as(:), cs(:)

    rhs = self%bulk_stencil%apply_along(f)
  end function diff_periodic

  pure function reverse(x)
    real, intent(in) :: x(:)
    real, allocatable :: reverse(:)
    reverse = x(size(x):1:-1)
  end function reverse

end module differentiate
