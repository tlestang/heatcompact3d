program main_field
  use field_module, only: Field
  implicit none

  real :: u0(16,16,16)

  type(Field) :: u
  u0 = 0.
  u = Field(u0)
end program main_field
