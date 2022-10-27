module coefficients_module
  real, parameter :: dirichlet_coeffs(8) = [ &
       & -5. / 2., 2., 0.5, 0., &
       & -3. / 4., 0., 3. / 4., 0. &
       & ]
  real, parameter :: neumann_odd_coeffs(8) = [ &
       & 0., 0., 0., 0., &
       & -7. / 9., -1. / 36., 7. / 9., 1. / 36. &
       & ]
  real, parameter :: sixth_order_compact_coeffs(4) = [ &
       - 1. / 36., &
       & - 7. / 9., &
       & + 7. / 9., &
       & + 1. / 36. &
       & ]
  real, parameter :: alpha_coeffs = 1. / 3.

  real, parameter :: dirichlet_2_coeffs(8) = [ &
       & 13., -27., 15., -1., &
       & 6. / 5., -12. / 5., 6. / 5., 0. &
       & ]
  real, parameter :: sixth_order_compact_2_coeffs(5) = [ &
       & 3. / 44., &
       & 12 / 11., &
       & -2. * (12. / 11. + 3. / 44.), &
       & 12. / 11., &
       & 3. / 44. &
       & ]
  real, parameter :: alpha_2_coeffs = 2. / 11.
  real, parameter :: gamma_coeffs = -1.
end module coefficients_module
