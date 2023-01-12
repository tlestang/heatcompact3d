Prototype for 3D heat equation solver

```math
\frac{\partial T}{\partial t}(x, y ,z) = \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} + \frac{\partial^2 T}{\partial z^2}
```
## Example

```fortran
use field, only: field_type
use time_integration: AB2_integrator_type

temp_field = field_type(u0, dx)
AB2 = AB2_integrator_type(starttime=0., endtime=0.3, dt)

call AB2%integrate(temp_field)

temp_field%dump("output.dat")
```

## Building and run the tests

Configure with

```
cmake -S . -B build
```

Tests are built and run from the build directory:

```
make && make test
```

## Generating and viewing the documentation

The documentation is generated using
[ford](https://forddocs.readthedocs.io/en/latest/). It is deployed
automatically from the `main` branch at
[https://tlestang.github.io/heatcompact3d/](tlestang.github.io/heatcompact3d/).

To preview the docs locally, you'll need `ford`:

```bash
pip install --user ford
```

Then, from the root of the repository:

```bash
ford heatcompact3d.md
```

You can then view the documentation in your browser, for instance:

```bash
firefox doc/index.html
```

## Layout

### High-level structure

Two client facing derived-types: `integrator_type` and `field_type`.

Derived type `field_type` (`module_field.f90`) models a three-dimensional
temperature field satisfying the heat equation discretized over a
unifomr spatial mesh of spacing `dx`. Type-bound function `rhs` of
`field_type` returns a new `field_type` instance evaluating the
right-hand side of the heat equation ($\Delta T$) over the mesh.

Type `time_integrator_type` (`module_time_integrator.f90`) provides a base abstract
derived type for time integration methods such as Euler,
Adams-Bashforth or Runge-Kutta methods. Type-bound method `integrate`
of `time_integrator_type` takes a `field_type` instance and overwrites it
with the result of the time integration. Currently available integrators:

- `euler_integrator_type`
- `AB2_integrator_type`
- `RK3_integrator_type`

### Spatial differentiation

Spatial differentation for both first order and second order
derivatives is implemented in the `differentiate` module
(`module_differentiate.f90`). The `differentiate` module provides access
to derived type `differentiator_type`, which `diff` type-bound
function implements differentiation of a one dimensional pencil, /e.g/

```math
g(i) = af_{i-2} + bf_{i-1} + af_{i+1} + bf_{i+2}
```

Instances of `differentiator_type` are created through constructors
`<method>_<order>`. Currently available differentiatiors are:

- `sixth_order_compact_{1,2}`

#### Differentiation stencils

A differentiation stencil of size $m$ takes the form

```math
ug_{i+1} + g + lg_{i-1} = \sum_{k=1}^{m} w_k f_{i+s_k}
```

it is modeled by derived-type `stencil_type` with the following components

- `nodes` (${s_k}_{1 \ leq k \leq m}$) :: stencil node's position
  relative to point to which the stencil is applied.
- `coeffs` (${w_k}_{1 \ leq k \leq m}$)
- `upper`, `lower` ($u$, $l$) coefficients forming the lower and upper
  diagonals of the tridiagnonal matrix describing the system.
  
A instance of `stencil_type` can be applied along a one dimensional
pencil

```fortran
s = stencil_type(...)
f, rhs = real(*)
rhs = s%apply_along(f)
```

#### Boundary schemes

Module `boundary` (`boundary_mod.f90`) composes four instances of
`stencil_type`.

- East boundary stencil (first order derivative)
- West boundary stencil (first order derivative)
- East boundary stencil (second order derivative)
- West boundary stencil (second order derivative)

The `boundary` provides access to specific boundary conditions through
special-purpose constructors. Currently implemented boundary schemes are

- `get_dirichlet_boundary` :: Returns a `boundary_type` instance
  describing a dirichlet boundary condition.

