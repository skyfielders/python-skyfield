"""Accuracy tests against data pulled from HORIZONS."""

from numpy import max
from skyfield import api
from skyfield.constants import AU_M
from skyfield.io import download
from skyfield.jpllib import Kernel

one_second = 1.0 / 24.0 / 60.0 / 60.0
arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
meter = 1.0 / AU_M

base = 'http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk'
de430_url = base + '/planets/de430.bsp'
de431_url = base + '/planets/de431.bsp'
jup310_url = base + '/satellites/jup310.bsp'

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
    k = Kernel(download(jup310_url))
    a = k.earth.geometry_of(k.callisto).at(tdb=2471184.5)
    compare(a.position.au,
      [-4.884815926454119E+00, -3.705745549073268E+00, -1.493487818022234E+00],
      0.001 * meter)
    compare(a.velocity.au_per_d,
      [9.604665478763035E-03, -1.552997751083403E-02, -6.678445860769302E-03],
      0.000001 * meter)

def test_callisto_astrometric():
    k = Kernel(download(jup310_url))
    a = k.earth.observe(k.callisto).at(utc=(2053, 10, 9))
    ra, dec, distance = a.radec()
    compare(ra._degrees, 217.1839292, 0.001 * arcsecond)
    compare(dec.degrees, -13.6892791, 0.001 * arcsecond)
    compare(distance.au, 6.31079291776184, 0.1 * meter)

def test_boston_geometry():
    k = Kernel(download(de430_url))
    jd = api.JulianDate(tdb=(2015, 3, 2), delta_t=67.185390 + 0.5285957)
    boston = api.Topos((42, 21, 24.1), (-71, 3, 24.8), x=0.003483, y=0.358609)
    a = boston.geometry_of(k.earth).at(jd)
    compare(a.position.km,
      [-1.764697476371664E+02, -4.717131288041386E+03, -4.274926422016179E+03],
      0.0027)  # TODO: try to get this < 1 meter

def test_moon_from_boston_geometry():
    k = Kernel(download(de430_url))
    jd = api.JulianDate(tdb=(2015, 3, 2), delta_t=67.185390 + 0.5285957)
    boston = api.Topos((42, 21, 24.1), (-71, 3, 24.8), x=0.003483, y=0.358609)
    a = boston.geometry_of(k.moon).at(jd)
    compare(a.position.au,
      [-1.341501206552443E-03, 2.190483327459023E-03, 6.839177007993498E-04],
      1.7 * meter)  # TODO: improve this

def test_moon_from_boston_astrometric():
    k = Kernel(download(de430_url))
    jd = api.JulianDate(tdb=(2015, 3, 2), delta_t=67.185390 + 0.5285957)
    boston = api.Topos((42, 21, 24.1), (-71, 3, 24.8), x=0.003483, y=0.358609)
    a = boston.observe(k.moon).at(jd)
    ra, dec, distance = a.radec()
    compare(ra._degrees, 121.4796470, 0.001 * arcsecond)
    compare(dec.degrees, 14.9108450, 0.001 * arcsecond)
    compare(distance.au, 0.00265828588792, 1.4 * meter)  # TODO: improve this
