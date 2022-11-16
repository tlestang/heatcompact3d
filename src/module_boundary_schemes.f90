module boundary_schemes
  use stencil, only: stencil_type
  implicit none

  type boundary_type
     type(stencil_type) :: first_order_east(2), second_order_east(2)
     type(stencil_type) :: first_order_west(2), second_order_west(2)
  end type boundary_type

  real, parameter :: sixth_order_compact_coeffs(4) = [ &
       & - 1. / 36., &
       & - 7. / 9., &
       & + 7. / 9., &
       & + 1. / 36. &
       & ]

  real, parameter :: sixth_order_compact_coeffs_2(5) = [ &
       & 3. / 44., &
       & 12. / 11., &
       & -2. * (12. / 11. + 3. / 44.), &
       & 12. / 11., &
       & 3. / 44. &
       & ]

contains

  pure type(stencil_type) function sixth_order_compact_stencil()
    sixth_order_compact_stencil = stencil_type( &
         & nodes = [-2, -1, 1, 2], &
         & coeffs = sixth_order_compact_coeffs, &
         & upper = 1. / 3.,  &
         & lower = 1. / 3. &
         & )
  end function sixth_order_compact_stencil

  pure type(stencil_type) function sixth_order_compact_second_stencil()
    sixth_order_compact_second_stencil = stencil_type( &
         & nodes = [-2, -1, 0, 1, 2], &
         & coeffs = sixth_order_compact_coeffs_2, &
         & upper = 2. / 11., &
         & lower = 2. / 11. &
         & )
  end function sixth_order_compact_second_stencil

  type(boundary_type) function get_dirichlet_boundary()
    type(stencil_type) :: first_order_east(2), second_order_east(2)
    type(stencil_type) :: first_order_west(2), second_order_west(2)

    first_order_east(1) = stencil_type( &
         & nodes = [0, 1, 2, 3], &
         & coeffs = [-5. / 2., 2., 0.5, 0.], &
         & lower = 0., upper = 2. &
         & )
    first_order_east(2) = stencil_type( &
         & nodes = [-1, 0, 1, 2], &
         & coeffs = [-3. / 4., 0., 3. / 4., 0.], &
         & lower = 1. / 4., upper = 1. / 4. &
         & )
    second_order_east(1) = stencil_type( &
         & nodes = [0, 1, 2, 3], &
         & coeffs = [13., -27., 15., -1.], &
         & lower = 0., upper = 11. &
         & )
    second_order_east(2) = stencil_type( &
         & nodes = [-1, 0, 1, 2], &
         & coeffs = [6. / 5., -12. / 5., 6. / 5., 0.], &
         & lower = 1. / 10., upper = 1. / 10. &
         & )
    first_order_west = first_order_east%flip() * (-1.)
    second_order_west = second_order_east%flip()
    get_dirichlet_boundary = boundary_type( &
         & first_order_east, second_order_east, first_order_west, second_order_west)
  end function get_dirichlet_boundary
end module boundary_schemes
