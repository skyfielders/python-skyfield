"""Accuracy tests against data pulled from HORIZONS."""

from numpy import max
from skyfield import api
from skyfield.constants import AU_M
from skyfield.jpllib import Kernel

one_second = 1.0 / 24.0 / 60.0 / 60.0
arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
meter = 1.0 / AU_M

def compare(value, expected_value, epsilon):
    if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
        assert max(abs(value - expected_value)) <= epsilon
    else:
        assert abs(value - expected_value) <= epsilon

def test_jupiter1():
    astrometric = api.sun(utc=(1980, 1, 1, 0, 0)).observe(api.jupiter)
    hlat, hlon, d = astrometric.ecliptic_latlon()
    compare(hlat.degrees, 1.013, 0.001)
    compare(hlon.degrees, 151.3229, 0.001)

def test_callisto_geometry():
    k = Kernel(open('jup310.bsp', 'rb'))
    a = k.earth.geometry_of(k.callisto).at(tdb=2471184.5)
    compare(a.position.au,
      [-4.884815926454119E+00, -3.705745549073268E+00, -1.493487818022234E+00],
      0.001 * meter)
    compare(a.velocity.au_per_d,
      [9.604665478763035E-03, -1.552997751083403E-02, -6.678445860769302E-03],
      0.000001 * meter)

def test_callisto_astrometric():
    k = Kernel(open('jup310.bsp', 'rb'))
    a = k.earth.observe(k.callisto).at(utc=(2053, 10, 9))
    ra, dec, distance = a.radec()
    print(ra._degrees)
    print(217.1839292)
    compare(ra._degrees, 217.1839292, 0.001 * arcsecond)
    compare(dec.degrees, -13.6892791, 0.001 * arcsecond)
    compare(distance.au, 6.31079291776184, 0.1 * meter)
