from numpy import abs

from skyfield.api import load
from skyfield.toposlib import Topos

angle = (15, 25, 35, 45)

def ts():
    yield load.timescale()

def test_beneath(ts, angle):
    t = ts.utc(2018, 1, 19, 14, 37, 55)
    # An elevation of 0 is more difficult for the routine's accuracy
    # than a very large elevation.
    top = Topos(latitude_degrees=angle, longitude_degrees=0, elevation_m=0)
    p = top.at(t)
    b = p.subpoint()
    error_degrees = abs(b.latitude.degrees - angle)
    error_mas = 60.0 * 60.0 * 1000.0 * error_degrees
    #print(b.latitude.degrees, deg, error_mas)
    assert error_mas < 0.1
