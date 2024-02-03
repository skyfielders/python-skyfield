from skyfield.api import load, wgs84
from numpy import array, sqrt

# Now let's try the SPICE 'surfpt.f' approach, and actually find the
# point of interception, while treating the Earth as an ellipsoid.
# Here are two helper routines:

# def VPERP(a, b):
#     return a - VPROJ(a, b)

from skyfield.functions import dots, length_of

def VPROJ(a, b):
    return dots(a,b) / dots(b,b) * b

assert list(VPROJ(array([6.0,  6.0, 6.0]),
             array([2.0,  0.0,  0.0]))) == ([6.0,  0.0,  0.0])

assert list(VPROJ(array([6.0,  6.0, 6.0]),
             array([-3.0,  0.0,  0.0]))) == ([6.0, -0.0, -0.0])

assert list(VPROJ(array([6.0,  6.0, 0.0]),
             array([0.0,  7.0,  0.0]))) == ([0.0,  6.0,  0.0])

assert list(VPROJ(array([6.0,  0.0, 0.0]),
             array([0.0,  0.0,  9.0]))) == ([0.0,  0.0,  0.0])

# And now, start work with the real position.

def main():
    ts = load.timescale()
    eph = load('de421.bsp')
    sun = eph['sun']
    moon = eph['Moon']
    earth = eph['earth']

    t = ts.utc(2017, 8, 21, 18, 30)

    e = earth.at(t)
    s = e.observe(sun).apparent()
    m = e.observe(moon).apparent()

    A = B = wgs84.radius.km  # x and y axes point to equator 6378.0
    C = wgs84.polar_radius.km

    print('Intersection, in ITRF coordinates:')
    shadow_itrf = line_and_ellipsoid_intersection(
        s.frame_xyz(itrs).km,
        m.frame_xyz(itrs).km - s.frame_xyz(itrs).km,
        array([A, B, C]),
    )
    print(shadow_itrf)

    # So where's that on Earth?
    from skyfield.positionlib import ICRF, Distance, Velocity
    position = ICRF.from_time_and_frame_vectors(t, itrs, Distance.km(shadow_itrf),
                                                Velocity(0))
    position.center = 399
    g = wgs84.geographic_position_of(position)
    print(g)

from skyfield.framelib import itrs

def line_and_ellipsoid_intersection(POSITN, U, radii):
    # And we're off to the races!  First, switch to unit-sphere mode.
    X = U / radii
    Y = POSITN / radii

    # Component of Y that's orthogonal to X?  Hopefully this all makes sense
    # later when I read back through it.
    P = Y - VPROJ( Y, X )

    # 'Find the magnitudes of Y and P.'
    YMAG = length_of(Y)
    PMAG = length_of(P)

    if YMAG > 1.0:  # y outside sphere
        if PMAG > 1.0:  # p also outside sphere
            exit('no intersection')
        YPROJ = Y - P
        if dots(YPROJ, X) > 0.0:
            exit('ray points away from sphere')
        if PMAG == 1.0:
            return P * radii
        # Must be a non-trivial intersection.
        sign = -1.0
    elif YMAG == 1.0:
        return POSITN
    else:  # y is inside sphere!
        sign = 1.0

    SCALE = sqrt(max(0.0, 1.0 - PMAG*PMAG))
    UX = X / length_of(X)
    POINT = P + sign * SCALE * UX
    return POINT * radii

main()

# Output:
# 'WGS84 latitude +36.2101 N longitude -85.9654 E elevation -3.1 m'
# Boom :)
