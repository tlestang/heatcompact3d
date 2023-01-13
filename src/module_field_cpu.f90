module field_cpu
  use field, only: field_type
  implicit none

  type, extends(field_type) :: field_cpu_type
   contains
     procedure, public :: rhs
  end type field_cpu_type

  interface field_cpu_type
     module procedure field_cpu_constructor
  end interface field_cpu_type

contains

  function field_cpu_constructor(initial, dx) result(afield)
    real, intent(in) :: initial(:, :, :) !! Initial state
    real, intent(in) :: dx !! Spatial mesh spacing
    type(field_cpu_type) :: afield
    allocate(afield%data, source=initial)
    afield%dx = dx
  end function field_cpu_constructor

  pure function rhs(self)
    !! Evaluates right hand side of heat equation on a field instance.
    !! \[ F(T) = \Delta T = \frac{\partial^2 T}{\partial x^2} +
    !! \frac{\partial^2 T}{\partial y^2} + \frac{\partial^2
    !! T}{\partial z^2} \]
    !!
    use differentiate, only: differentiator_type, &
         & sixth_order_compact_1, sixth_order_compact_2

    class(field_cpu_type), intent(in) :: self
    real, allocatable :: ddx(:, :, :), ddy(:, :, :), ddz(:, :, :)
    class(field_type), allocatable :: rhs
    integer :: ix, iy, iz
    real :: dx2
    class(differentiator_type), allocatable :: differ

    dx2 = self%dx * self%dx
    differ = sixth_order_compact_2() ! Periodic boundaries

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

    allocate(rhs, source=self)
    rhs%data = ddx + ddy + ddz
    rhs%dx = self%dx
  end function rhs

end module field_cpu
