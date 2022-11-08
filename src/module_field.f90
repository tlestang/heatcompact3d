module field_module
  implicit none

  type :: Field
     !! Implement a 3D scalar field, for instance a temperature field.
     !! ```f90
     !! type(field_type) = afield
     !! real :: u0(nx, ny, nz), dx
     !! afield = field_type(u0, dx)
     !! ```
     private
     real, allocatable :: data(:, :, :)
     real :: dx
     !! Discrete mesh spacing
   contains
     procedure, public :: nx, ny, nz
     procedure, public :: is_equal
     procedure, public :: rhs
     procedure, private :: field_add_field, field_sub_field, &
          & field_mul_real
     generic :: operator(+) => field_add_field
     generic :: operator(-) => field_sub_field
     generic :: operator(*) => field_mul_real
  end type Field

  interface Field
     module procedure constructor
  end interface Field

contains

  function constructor(initial, dx) result(afield)
    real, intent(in) :: initial(:, :, :) !! Initial state
    real, intent(in) :: dx !! Spatial mesh spacing
    type(Field) :: afield
    allocate(afield%data, source=initial)
    afield%dx = dx
  end function constructor

  pure integer function nx(self)
    !! Returns domain size in \(\mathbf{x}\) direction
    class(Field), intent(in) :: self
    nx = size(self%data, 1)
  end function nx

  pure integer function ny(self)
    !! Returns domain size in \(\mathbf{y}\) direction
    class(Field), intent(in) :: self
    ny = size(self%data, 2)
  end function ny

  pure integer function nz(self)
    !! Returns domain size in \(\mathbf{z}\) direction
    class(Field), intent(in) :: self
    nz = size(self%data, 3)
  end function nz

  pure logical function is_equal(self, lhs, tol)
    !! Compare two field_type instance based on their data value
    !! ```f90
    !! f1 = field_type(u0, dx)
    !! f2 = field_type(u0, dx2)
    !! f1%is_equal(f2) ! true
    !! ```
    class(Field), intent(in) :: self !! Right hand side of comparison
    class(Field), intent(in) :: lhs !! Left hand side of comparison
    real, intent(in) :: tol !! Absolute tolerance when comparing
                            !! fields values
    logical, allocatable :: elmt_is_equal(:, :, :)

    elmt_is_equal = abs(self%data - lhs%data) < tol
    is_equal = all(elmt_is_equal)
  end function is_equal

  pure function rhs(self)
    !! Evaluates right hand side of heat equation on a field instance.
    !! \[ F(T) = \Delta T = \frac{\partial^2 T}{\partial x^2} +
    !! \frac{\partial^2 T}{\partial y^2} + \frac{\partial^2
    !! T}{\partial z^2} \]
    !!
    use differentiate, only: differentiator_type, &
       & sixth_order_compact, sixth_order_compact_2nd

    class(Field), intent(in) :: self
    real, allocatable :: ddx(:, :, :), ddy(:, :, :), ddz(:, :, :)
    type(Field) :: rhs
    integer :: ix, iy, iz
    real :: dx2
    class(differentiator_type), allocatable :: differ

    dx2 = self%dx * self%dx
    differ = sixth_order_compact_2nd()

    allocate(ddx, source=self%data)
    do iz = 1,self%nz()
       do iy = 1,self%ny()
          ddx(:, iy, iz) = differ%diff(self%data(:, iy, iz), dx2)
       end do
    end do

    allocate(ddy, source=self%data)
    do iz = 1,self%nz()
       do ix = 1,self%nx()
          ddy(ix, :, iz) = differ%diff(self%data(ix, :, iz), dx2)
       end do
    end do

    allocate(ddz, source=self%data)
    do iy = 1,self%ny()
       do ix = 1,self%nx()
          ddz(ix, iy, :) = differ%diff(self%data(ix, iy, :), dx2)
       end do
    end do

    rhs%data = ddx + ddy + ddz
    rhs%dx = self%dx
  end function rhs

  pure type(Field) function field_add_field(self, afield)
    class(Field), intent(in) :: self, afield
    field_add_field%data = self%data + afield%data
    field_add_field%dx = self%dx
  end function field_add_field

  pure type(Field) function field_sub_field(self, afield)
    class(Field), intent(in) :: self, afield
    field_sub_field%data = self%data - afield%data
    field_sub_field%dx = self%dx
  end function field_sub_field

  pure type(Field) function field_mul_real(self, a)
    !! Multiply a `field_type` instance by a `real` number.
    !! ```
    !! f1 = field_type(u0, dx)
    !! f2 = f1 * 1.3
    !! f2%is_equal(field_type(u0 * 1.3, dx)) ! true
    !! ```
    class(Field), intent(in) :: self !! Left hand side
    real, intent(in) :: a !! Scalar to multiply field instance with
    field_mul_real%data = self%data * a
    field_mul_real%dx = self%dx
  end function field_mul_real

end module field_module
