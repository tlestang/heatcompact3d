module field_gpu
  use field, only: field_type
  implicit none

  type, extends(field_type) :: field_gpu_type
     real, allocatable, device :: data(:, :, :)
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

    allocate(rhs, mold=self)
    rhs%data = 2. * self%data
    rhs%dx = self%dx
  end function rhs
end module field_gpu
