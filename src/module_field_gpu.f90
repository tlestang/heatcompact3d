module field_gpu
  use cudafor
  use field, only: field_type
  implicit none

  type, extends(field_type) :: field_gpu_type
     real, allocatable, device :: data_dev(:, :, :)
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
    class(field_type), allocatable :: rhs
    type(dim3) :: threads

    threads = dim3(16, 16, 16)
    allocate(field_gpu_type :: rhs)
    select type (rhs)
       type is (field_gpu_type)
          allocate(rhs%data_dev, mold=self%data_dev)
          call rhs_kernel<<<1, threads>>>(self%data_dev, rhs%data_dev)
       end select
     end function rhs

  attributes(global) pure subroutine rhs_kernel(u, rhs)
    real, intent(in) :: u(:, :, :)
    real, intent(out) :: rhs(:, :, :)
    integer :: i, j, k

    i = threadIdx%x
    j = threadIdx%y
    k = threadIdx%z
    rhs(i, j, k) = 2. * u(i, j, k)
   end subroutine rhs_kernel

end module field_gpu
