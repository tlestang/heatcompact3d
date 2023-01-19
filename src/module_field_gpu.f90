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
       ! x-directional derivatives
       grid = dim3(ceiling(real(self%ny())/th), ceiling(real(self%nz())/th), 1)
       threads = dim3(th, th, 1)
       call rhs_x<<<grid, threads>>>(data_dev, rhs_dev, temp_dev)

       ! y-directional derivatives
       grid = dim3(ceiling(real(self%nx())/th), ceiling(real(self%nz())/th), 1)
       threads = dim3(th, th, 1)
       call rhs_y<<<grid, threads>>>(data_dev, rhs_dev, temp_dev)

       ! z-directional derivatives
       grid = dim3(ceiling(real(self%nx())/th), ceiling(real(self%ny())/th), 1)
       threads = dim3(th, th, 1)
       call rhs_z<<<grid, threads>>>(data_dev, rhs_dev, temp_dev)

       rhs = rhs_dev
       end select
     end function rhs

  attributes(global) pure subroutine rhs_x(data, rhs, temp)
    real, intent(in) :: data(:, :, :)
    real, intent(out) :: rhs(:, :, :)
    real, intent(inout) :: temp(:, :, :)
    integer :: x, y, z, nx, ny, nz
    real :: dx

    nx = size(data, dim=1)
    ny = size(data, dim=2)
    nz = size(data, dim=3)
    y = blockDim%x*(blockIdx%x-1) + threadIdx%x
    z = blockDim%y*(blockIdx%y-1) + threadIdx%y

    dx = 2. * acos(-1.) / (nx - 1)

    if ( y <= ny .and. z <= nz ) then
       ! Construct the RHS and apply Forward Thomas
       do x = 1, nx
          temp(x, y, z) = data(x, y, z)
       end do
       ! Backward Thomas
       do x = nx, 1, -1
          temp(x, y, z) = -1. * cos((x-1)*dx) * cos((y-1)*dx) * cos((z-1)*dx) !2.*temp(x, y, z)
       end do
       ! Periodic last pass
       do x = 1, nx
          rhs(x, y, z) = temp(x, y, z)
       end do
    end if

   end subroutine rhs_x

  attributes(global) pure subroutine rhs_y(data, rhs, temp)
    real, intent(in) :: data(:, :, :)
    real, intent(inout) :: rhs(:, :, :)
    real, intent(inout) :: temp(:, :, :)
    integer :: x, y, z, nx, ny, nz
    real :: dx

    nx = size(data, dim=1)
    ny = size(data, dim=2)
    nz = size(data, dim=3)
    x = blockDim%x*(blockIdx%x-1) + threadIdx%x
    z = blockDim%y*(blockIdx%y-1) + threadIdx%y

    dx = 2. * acos(-1.) / (nx - 1)

    if ( x <= nx .and. z <= nz ) then
       ! Construct the RHS and apply Forward Thomas
       do y = 1, ny
          temp(x, y, z) = data(x, y, z)
       end do
       ! Backward Thomas
       do y = ny, 1, -1
          temp(x, y, z) = -1. * cos((x-1)*dx) * cos((y-1)*dx) * cos((z-1)*dx) !2.*temp(x, y, z)
       end do
       ! Periodic last pass
       do y = 1, ny
          rhs(x, y, z) = rhs(x, y, z) + temp(x, y, z)
       end do
    end if

   end subroutine rhs_y

  attributes(global) pure subroutine rhs_z(data, rhs, temp)
    real, intent(in) :: data(:, :, :)
    real, intent(inout) :: rhs(:, :, :)
    real, intent(inout) :: temp(:, :, :)
    integer :: x, y, z, nx, ny, nz
    real :: dx

    nx = size(data, dim=1)
    ny = size(data, dim=2)
    nz = size(data, dim=3)
    x = blockDim%x*(blockIdx%x-1) + threadIdx%x
    y = blockDim%y*(blockIdx%y-1) + threadIdx%y

    dx = 2. * acos(-1.) / (nx - 1)

    if ( x <= nx .and. y <= ny ) then
       ! Construct the RHS and apply Forward Thomas
       do z = 1, nz
          temp(x, y, z) = data(x, y, z)
       end do
       ! Backward Thomas
       do z = nz, 1, -1
          temp(x, y, z) = -1. * cos((x-1)*dx) * cos((y-1)*dx) * cos((z-1)*dx) !2.*temp(x, y, z)
       end do
       ! Periodic last pass
       do z = 1, nz
          rhs(x, y, z) = rhs(x, y, z) + temp(x, y, z)
       end do
    end if

   end subroutine rhs_z
  
end module field_gpu
