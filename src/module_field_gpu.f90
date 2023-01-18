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
    real, allocatable, device :: rhs_dev(:, :, :), data_dev(:, :, :), temp_dev(:, :, :)
    type(dim3) :: threads, grid
    integer, parameter :: th = 16

    select type (self)
    type is (field_gpu_type)
       allocate(data_dev, mold=self%data)
       allocate(rhs_dev, mold=self%data)
       allocate(temp_dev, mold=self%data)

       rhs_dev = self%data
       grid = dim3(ceiling(real(self%ny())/th), ceiling(real(self%nz())/th), 1)
       threads = dim3(th, th, 1)
       call rhs_x<<<grid, threads>>>(data_dev, rhs_dev, temp_dev)

       rhs = rhs_dev
       end select
     end function rhs

  attributes(global) pure subroutine rhs_x(data, rhs, temp)
    real, intent(in) :: data(:, :, :)
    real, intent(out) :: rhs(:, :, :)
    real, intent(inout) :: temp(:, :, :)
    integer :: x, y, z, nx, ny, nz

    nx = size(data, dim=1)
    ny = size(data, dim=2)
    nz = size(data, dim=3)
    y = blockDim%x*(blockIdx%x-1) + threadIdx%x
    z = blockDim%y*(blockIdx%y-1) + threadIdx%y

    if ( y <= ny .and. z <= nz ) then
       do x = 1, nx
          temp(x, y, z) = data(x, y, z)
       end do
       do x = 1, nx
          temp(x, y, z) = 2.*temp(x, y, z)
       end do
       do x = 1, nx
          rhs(x, y, z) = temp(x, y, z)
       end do
    end if

   end subroutine rhs_x
  
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
