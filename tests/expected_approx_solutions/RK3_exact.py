import sympy
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.vector import CoordSys3D
from sympy.vector import divergence as div
from sympy.vector import gradient as grad
from sympy import Rational

h = sympy.symbols('h')
C = CoordSys3D('C')
u0  = cos(C.x) * cos(C.y) * cos(C.z)

for i in [1, 2]:
    u1 = u0 + Rational('8/15') * h * div(grad(u0))
    u2 = u1 + Rational('5/12') * h * div(grad(u1)) - Rational('17/60') * h * div(grad(u0))
    u3 = u2 + Rational('3/4') * h * div(grad(u2)) - Rational('5/12') * h * div(grad(u1))

    u0 = u3
print(sympy.factor(u3))

