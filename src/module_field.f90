module field_module
  implicit none

  type :: Field
     private
     real, allocatable :: data(:, :, :)
   ! contains
   !   procedure, public :: rhs
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

  ! pure function rhs(self)
  !   type(Field), intent(in) :: self
  !   real, allocatable :: d2_dx2, d2_dy2, d2_dz2
  !   integer :: nx, ny, nz
  !   type(Field) :: rhs

  !   nx = size(self%data, 1)
  !   ny = size(self%data, 2)
  !   nz = size(self%data, 3)

  !   allocate(d2_dx2, source=self%data)
  !   do iz = 1,nz
  !      do iy = 1,nz
  !         d2_dx2(:, iy, iz) = diff2(data(:, iy, iz))
  !      end do
  !   end do

  !   allocate(d2_dy2, source=self%data)
  !   do iz = 1,nz
  !      do ix = 1,nx
  !         d2_dy2(ix, :, iz) = diff2(data(ix, :, iz))
  !      end do
  !   end do

  !   allocate(d2_dz2, source=self%data)
  !   do iy = 1,ny
  !      do ix = 1,nx
  !         d2_dz2(ix, iy, :) = diff2(data(ix, iy, :))
  !      end do
  !   end do

  !   rhs%data = d2_dx2 + d2_dy2 + d2_dz2
  ! end function rhs
end module field_module
