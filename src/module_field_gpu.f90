module field_gpu
  use cudafor
  use field, only: field_type
  implicit none

  type, extends(field_type) :: field_gpu_type
   contains
     procedure, public :: rhs
  end type field_gpu_type

  interface field_gpu_type
     module procedure field_gpu_constructor
  end interface field_gpu_type

contains

  function field_gpu_constructor(initial, dx) result(afield)
    real, intent(in) :: initial(:, :, :) !! Initial state
    real, intent(in) :: dx !! Spatial mesh spacing
    type(field_gpu_type) :: afield
    allocate(afield%data, source=initial)
    afield%dx = dx
  end function field_gpu_constructor

  pure function rhs(self)
    class(field_gpu_type), intent(in) :: self
    real, allocatable :: rhs(:, :, :)
    real, allocatable, device :: rhs_dev(:, :, :), data_dev(:, :, :)
    type(dim3) :: threads

    threads = dim3(16, 16, 16)
    select type (self)
    type is (field_gpu_type)
       allocate(data_dev, mold=self%data)
       allocate(rhs_dev, mold=self%data)
       rhs_dev = self%data
       call rhs_kernel<<<1, threads>>>(data_dev, rhs_dev)
       rhs = rhs_dev
       end select
     end function rhs

  attributes(global) pure subroutine rhs_kernel(data, rhs)
    real, intent(in) :: data(:, :, :)
    real, intent(out) :: rhs(:, :, :)
    integer :: i, j, k

    i = threadIdx%x
    j = threadIdx%y
    k = threadIdx%z
    rhs(i, j, k) = 2. * data(i, j, k)
   end subroutine rhs_kernel

end module field_gpu
