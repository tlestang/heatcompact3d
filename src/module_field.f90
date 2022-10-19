module field_module
  implicit none

  type :: Field
     private
     real, allocatable :: data(:, :, :)
   contains
     procedure, public :: nx, ny, nz
     procedure, public :: is_equal
     procedure, public :: rhs
  end type Field

  interface Field
     module procedure constructor
  end interface Field

contains

  function constructor(initial) result(afield)
    real, intent(in) :: initial(:, :, :)
    type(Field) :: afield
    allocate(afield%data, source=initial)
  end function constructor

  pure integer function nx(self)
    class(Field), intent(in) :: self
    nx = size(self%data, 1)
  end function nx

  pure integer function ny(self)
    class(Field), intent(in) :: self
    ny = size(self%data, 2)
  end function ny

  pure integer function nz(self)
    class(Field), intent(in) :: self
    nz = size(self%data, 3)
  end function nz

  pure logical function is_equal(self, lhs, tol)
    class(Field), intent(in) :: self
    class(Field), intent(in) :: lhs
    real, intent(in) :: tol
    logical, allocatable :: elmt_is_equal(:, :, :)

    elmt_is_equal = abs(self%data - lhs%data) < tol
    is_equal = all(elmt_is_equal)
  end function is_equal

  pure function rhs(self, dx)
    use differentiate, only: diff2

    real, intent(in) :: dx
    class(Field), intent(in) :: self
    real, allocatable :: ddx(:, :, :), ddy(:, :, :), ddz(:, :, :)
    type(Field) :: rhs
    integer :: ix, iy, iz

    allocate(ddx, source=self%data)
    do iz = 1,self%nz()
       do iy = 1,self%ny()
          ddx(:, iy, iz) = diff2(self%data(:, iy, iz), dx)
       end do
    end do

    allocate(ddy, source=self%data)
    do iz = 1,self%nz()
       do ix = 1,self%nx()
          ddy(ix, :, iz) = diff2(self%data(ix, :, iz), dx)
       end do
    end do

    allocate(ddz, source=self%data)
    do iy = 1,self%ny()
       do ix = 1,self%nx()
          ddz(ix, iy, :) = diff2(self%data(ix, iy, :), dx)
       end do
    end do

    rhs%data = ddx + ddy + ddz
  end function rhs
end module field_module
