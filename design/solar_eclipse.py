from skyfield.api import load, wgs84
from numpy import array, sqrt

# Now let's try the SPICE 'surfpt.f' approach, and actually find the
# point of interception, while treating the Earth as an ellipsoid.
# Here are two helper routines:

# def VPERP(a, b):
#     return a - VPROJ(a, b)

from skyfield.functions import dots, length_of

def vector_projection(a, b):
    return dots(a,b) / dots(b,b) * b

assert list(vector_projection(array([6.0,  6.0, 6.0]),
             array([2.0,  0.0,  0.0]))) == ([6.0,  0.0,  0.0])

assert list(vector_projection(array([6.0,  6.0, 6.0]),
             array([-3.0,  0.0,  0.0]))) == ([6.0, -0.0, -0.0])

assert list(vector_projection(array([6.0,  6.0, 0.0]),
             array([0.0,  7.0,  0.0]))) == ([0.0,  6.0,  0.0])

assert list(vector_projection(array([6.0,  0.0, 0.0]),
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

class A(object):
    __getitem__ = array
A = A()

import numpy as np
from numpy import nan #, sign

nan3 = A[nan,nan,nan]

def line_and_ellipsoid_intersection(line_start, line_direction, radii):
    # Based on `surfpt.f` from the SPICE Toolkit.
    # TODO: vectorize

    # Scale coordinates so the ellipsoid becomes the unit sphere.
    start = line_start / radii
    direction = line_direction / radii

    # What point on the line is closest to the sphere's center?
    closest_point = start - vector_projection(start, direction)

    startmag = length_of(start)
    pmag = length_of(closest_point)
    startPROJ = start - closest_point

    j = 2 + np.sign(startmag - 1.0).astype(int) * 2
    k = dots(startPROJ, direction) > 0.0
    i = j + k

    sign_table = A[
        +1,  # startmag < 1: looking out from inside sphere
        +1,  # startmag < 1: (same)
        -1,  # startmag = 1: on surface, so choose - to return `start`
        +1,  # startmag = 1: on surface, so choose + to return `start`
        -1,  # startmag > 1: outside sphere, looking towards it
        nan, # startmag > 1: outside sphere, but looking away from it
    ]
    # sign_table = A[
    #     +1,  # startmag < 1, same direction: looking out from inside sphere
    #     -1,  # startmag = 1, ...: on surface of sphere ...
    #     -1,  # startmag > 1, ...: outside sphere ...
    #     +1,  # startmag < 1, ...: inside sphere
    #     +1,  # startmag = 1, ...: on surface of sphere ...
    #     nan, # startmag > 1, ...: outside sphere looking away from it
    # ]
    sign = sign_table[i]

    # if startmag > 1.0:          # start is outside sphere
    #     startPROJ = start - closest_point
    #     if dots(startPROJ, direction) > 0.0:
    #         #return nan3
    #         sign = nan
    #     # if pmag > 1.0:
    #     #     return nan3
    #     # if pmag == 1.0:
    #     #     return closest_point * radii
    #     else:
    #         sign = -1.0
    # elif startmag == 1.0:       # start is on surface of sphere
    #     print('HERE')
    #     print('*******', startmag, pmag)
    #     # return line_start
    #     # sign = 1.0
    #     startPROJ = start - closest_point
    #     sign = np.sign(dots(startPROJ, direction))
    #     # if  > 0.0:
    #     #     sign = 1.0
    #     # else:
    #     #     sign = -1.0
    # else:                       # start is inside sphere
    #     sign = 1.0

    half_chord_length = sqrt(1.0 - pmag*pmag)  # TODO: issues warning?!
    unit_direction = direction / length_of(direction)
    print('*****',closest_point, sign, half_chord_length, unit_direction)
    intersection = closest_point + sign * half_chord_length * unit_direction
    return intersection * radii

def test():
    f = line_and_ellipsoid_intersection
    assert list(f(A[-5,0,0], A[1,0,0], A[1,2,3])) == [-1,0,0]
    assert list(f(A[-5,2,0], A[1,0,0], A[1,2,3])) == [0,2,0]
    assert list(f(A[5,0,3], A[-1,0,0], A[1,2,3])) == [0,0,3]

    starts = A[[-5,0,0], [-5,2,0], [0,-5,0], [0,-5,4], [0,0,-5], [1,0,-5]]
    surface_points = A[[1,0,0], [0,2,0], [0,0,4], [-1,0,0], [0,-2,0], [0,0,-4]]
    radii = A[1,2,4]

    def check(starts, answers):
        for start, answer in zip(starts, answers):
            actual = f(start * 1.0, direction * 1.0, radii)
            print(start, direction, answer, '?=', actual)
            assert str(answer * 1.0) == str(actual)

    direction = A[1,0,0]
    answers = A[[-1,0,0], [0,2,0], nan3, nan3, nan3, nan3]
    check(starts, answers)
    check(surface_points, surface_points)

    direction = A[0,1,0]
    answers = A[nan3, nan3, [0,-2,0], [0,0,4], nan3, nan3]
    check(starts, answers)
    check(surface_points, surface_points)

    direction = A[0,0,1]
    answers = A[nan3, nan3, nan3, nan3, [0,0,-4], [1,0,0]]
    check(starts, answers)
    check(surface_points, surface_points)

    direction = A[-1,0,0]
    answers = A[nan3, nan3, nan3, nan3, nan3, nan3]
    check(starts, answers)
    check(surface_points, surface_points)

test()
main()

# Output:
# 'WGS84 latitude +36.2101 N longitude -85.9654 E elevation -3.1 m'
# Boom :)
