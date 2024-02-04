from skyfield.api import load, wgs84
from numpy import array, nan

from skyfield.functions import dots
from skyfield.geometry import line_and_ellipsoid_intersection

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

    #t = ts.utc(2017, 8, 21, 18, 30)
    t = ts.utc(2017, 8, 21, 16, range(0, 300, 15))

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
    position = ICRF.from_time_and_frame_vectors(
        t, itrs, Distance.km(shadow_itrf), Velocity(0),
    )
    position.center = 399
    g = wgs84.geographic_position_of(position)
    print(g)

    latitude = g.latitude.degrees
    longitude = g.longitude.degrees
    midpoint = len(latitude) // 2

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    o = ccrs.Orthographic(
        central_longitude=longitude[midpoint],
        central_latitude=latitude[midpoint],
    )

    fig = plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=o)
    ax.coastlines(resolution='110m')
    ax.gridlines()
    ax.set_global()
    ax.plot(
        longitude, latitude,
        color='blue', linewidth=1, marker='+',
        transform=ccrs.Geodetic(),
    )

    fig.savefig('eclipse.png')

from skyfield.framelib import itrs

class A(object):
    __getitem__ = array
A = A()

nan3 = A[nan,nan,nan]

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
