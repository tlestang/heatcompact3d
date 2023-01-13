module field_gpu
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

  attributes(global) pure function rhs(self)
    class(field_gpu_type), intent(in) :: self
    class(field_type), allocatable :: rhs

    i = threadIdx%x
    j = threadIdx%y
    k = threadIdx%z
    allocate(rhs, mold=self)
    rhs%data_dev(i, j, k) = 2. * self%data_dev(i, j, k)
    rhs%dx = self%dx
  end function rhs

  pure subroutine dev_to_host(self)
    class(field_gpu_type), intent(inout) :: self
    self%data = self%data_dev
  end subroutine dev_to_host
end module field_gpu
