import pytest
from numpy import abs
from skyfield.api import earth, mars
from skyfield.jpllib import Ephemeris

try:
    import de405
    de405 = Ephemeris(de405)
except ImportError:
    pytestmark = pytest.mark.skipif(True, reason='de405 unavailable')

arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0

def compare(value, benchmark_value, tolerance):
    assert abs(value - benchmark_value) < tolerance

def test_mars0():
    a = de405.earth(tt=2440423.345833333).observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 16.0296606272219, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.12731030858147, 0.001 * arcsecond)

def test_mars1():
    a = de405.earth(tt=2448031.5).observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 23.545034875459514, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.882249043221036, 0.001 * arcsecond)

def test_mars2():
    a = de405.earth(tt=2451545.0).observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 22.034936616343344, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.180707411034978, 0.001 * arcsecond)

def test_mars3():
    a = de405.earth(tt=2456164.5).observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 13.894324196598355, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.122808318928705, 0.001 * arcsecond)

