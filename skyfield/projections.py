
from numpy import sqrt
from .functions import length_of

def _derive_stereographic():
    """Compute the formulae to cut-and-paste into the routine below."""
    from sympy import symbols, atan2, acos, rot_axis1, rot_axis3, Matrix
    x_c, y_c, z_c, x, y, z = symbols('x_c y_c z_c x y z')

    # The angles we'll need to rotate through.
    around_z = atan2(x_c, y_c)
    around_x = acos(-z_c)

    # Apply rotations to produce an "o" = output vector.
    v = Matrix([x, y, z])
    xo, yo, zo = rot_axis1(around_x) * rot_axis3(-around_z) * v

    # Which we then use the stereographic projection to produce the
    # final "p" = plotting coordinates.
    xp = xo / (1 - zo)
    yp = yo / (1 - zo)

    return xp, yp

def _optimize(expressions):
    from sympy import cse, numbered_symbols
    commons, outputs = cse(
        expressions,
        numbered_symbols('t'),
        optimizations='basic',
    )
    for symbol, expr in commons:
        print(symbol, '=', expr)
    print()
    for expr in outputs:
        print(expr)

def build_stereographic_projection(center):
    """Compute *x* and *y* coordinates at which to plot the positions."""

    # TODO: Computing the center should really be done using
    # optimization, as in:
    # https://math.stackexchange.com/questions/409217/
    p = center.xyz.au
    u = p / length_of(p)
    if len(u.shape) > 1:
        c = u.mean(axis=1)
        c = c / length_of(c)
    else:
        c = u
    x_c, y_c, z_c = c

    def project(position):
        p = position.xyz.au
        u = p / length_of(p)
        x, y, z = u
#        x_out = (x*y_c/sqrt(x_c**2 + y_c**2) - x_c*y/sqrt(x_c**2 + y_c**2))/(x*x_c*sqrt(-z_c**2 + 1)/sqrt(x_c**2 + y_c**2) + y*y_c*sqrt(-z_c**2 + 1)/sqrt(x_c**2 + y_c**2) + z*z_c + 1)
#        y_out = (-x*x_c*z_c/sqrt(x_c**2 + y_c**2) - y*y_c*z_c/sqrt(x_c**2 + y_c**2) + z*sqrt(-z_c**2 + 1))/(x*x_c*sqrt(-z_c**2 + 1)/sqrt(x_c**2 + y_c**2) + y*y_c*sqrt(-z_c**2 + 1)/sqrt(x_c**2 + y_c**2) + z*z_c + 1)
#        return x_out, y_out

        t0 = 1/sqrt(x_c**2 + y_c**2)
        t1 = x*x_c
        t2 = sqrt(-z_c**2 + 1)
        t3 = t0*t2
        t4 = y*y_c
        t5 = 1/(t1*t3 + t3*t4 + z*z_c + 1)
        t6 = t0*z_c

        return t0*t5*(x*y_c - x_c*y), -t5*(t1*t6 - t2*z + t4*t6)

    return project
