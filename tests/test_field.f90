program test_field
  use iso_fortran_env, only: stderr => error_unit
  use field_module, only: Field
  implicit none

  real :: u0(16, 16, 16)
  type(Field) :: temp_field, rhs, expected
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
  temp_field = Field(u0)

  allpass = .true.

  if (.not. all([nx, ny, nz] == &
       &[temp_field%nx(), temp_field%ny(), temp_field%nz()])) then
     write(stderr, '(a)') 'Field data has the correct rank and size... failed.'
     allpass = .false.
  else
     write(stderr, '(a)') 'Field data has the correct rank and size... passed.'
  end if

  if(.not. temp_field%is_equal(Field(u0), tol)) then
     write(stderr, '(a)') 'Field equal operator... failed.'
     allpass = .false.
  else
     write(stderr, '(a)') 'Field equal operator... passed.'
  end if

  expected = Field(-1. * 3. * u0)
  rhs = temp_field%rhs(dx)
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
