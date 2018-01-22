from numpy import abs

from skyfield.api import load
from skyfield.constants import AU_M, ERAD
from skyfield.positionlib import Geocentric
from skyfield.toposlib import Topos

def ts():
    yield load.timescale()

def test_beneath(ts):
    t = ts.utc(2018, 1, 19, 14, 37, 55)
    for deg in 15, 25, 35, 45:
        # An elevation of 0 is more difficult for the routine's accuracy
        # than a very large elevation.
        top = Topos(latitude_degrees=deg, longitude_degrees=0, elevation_m=0)
        p = top.at(t)
        b = Topos.subpoint_beneath(p)
        error_degrees = abs(b.latitude.degrees - deg)
        error_mas = 60.0 * 60.0 * 1000.0 * error_degrees
        #print(b.latitude.degrees, deg, error_mas)
        assert error_mas < 0.1
