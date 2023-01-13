program test_field_cuda
  use iso_fortran_env, only: stderr => error_unit
  use field_gpu, only: field_gpu_type
  implicit none

  real :: u0(16, 16, 16)
  class(field_gpu_type), allocatable :: temp_field, rhs, expected
  real, parameter :: tol = 0.1
  integer :: i, j, k, nx, ny, nz
  logical :: allpass
  real :: dx

  nx = size(u0, 1)
  ny = size(u0, 2)
  nz = size(u0, 3)
  dx = 2. * acos(-1.) / (nx - 1)

  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           u0(i, j, k) = &
                & cos((i-1)*dx) * cos((j-1)*dx) * &
                & cos((k-1)*dx)
        end do
     end do
  end do
  temp_field = field_gpu_type(u0, dx)

  allpass = .true.

  temp_field%data_dev = temp_field%data
  rhs = temp_field%rhs()
  rhs%data = rhs%data_dev

  expected = field_gpu_type(-1. * 3. * u0, dx)
  if(.not. expected%is_equal(rhs, tol)) then
     write(stderr, '(a)') 'Field right hand side is computed correctly... failed.'
     allpass = .false.
  else
     write(stderr, '(a)') 'Field right hand side is computed correctly... passed.'
  end if

  if (allpass) then
     write(stderr, '(a)') 'All tests passed successfully'
  else
     write(stderr, '(a)') '!!! Some tests failed !!!'
  end if
end program test_field
