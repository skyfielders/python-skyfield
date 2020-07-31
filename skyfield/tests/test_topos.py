from numpy import abs

from skyfield import constants
from skyfield.api import Topos, load
from skyfield.functions import length_of

angle = (-15, 15, 35, 45)

def ts():
    yield load.timescale()

def test_velocity():
    # It looks like this is a sweet spot for accuracy: presumably a
    # short enough fraction of a second that the vector does not time to
    # change direction much, but long enough that the direction does not
    # get lost down in the noise.
    factor = 300.0

    ts = load.timescale()
    t = ts.utc(2019, 11, 2, 3, 53, [0, 1.0 / factor])
    jacob = Topos(latitude_degrees=36.7138, longitude_degrees=-112.2169)
    p = jacob.at(t)
    velocity1 = p.position.km[:,1] - p.position.km[:,0]
    velocity2 = p.velocity.km_per_s[:,0]
    print(length_of(velocity2 - factor * velocity1))
    assert length_of(velocity2 - factor * velocity1) < 0.0007

def test_itrf_vector():
    ts = load.timescale()
    t = ts.utc(2019, 11, 2, 3, 53)
    top = Topos(latitude_degrees=45, longitude_degrees=0,
                elevation_m=constants.AU_M - constants.ERAD)
    x, y, z = top.at(t).itrf_xyz().au
    assert abs(x - 0.7071) < 1e-4
    assert abs(y - 0.0) < 1e-14
    assert abs(z - 0.7071) < 1e-4

def test_beneath(ts, angle):
    t = ts.utc(2018, 1, 19, 14, 37, 55)
    # An elevation of 0 is more difficult for the routine's accuracy
    # than a very large elevation.
    top = Topos(latitude_degrees=angle, longitude_degrees=angle, elevation_m=0)
    p = top.at(t)
    b = p.subpoint()

    error_degrees = abs(b.latitude.degrees - angle)
    error_mas = 60.0 * 60.0 * 1000.0 * error_degrees
    assert error_mas < 0.1

    error_degrees = abs(b.longitude.degrees - angle)
    error_mas = 60.0 * 60.0 * 1000.0 * error_degrees
    assert error_mas < 0.1
