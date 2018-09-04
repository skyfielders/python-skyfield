"""Accuracy tests against data pulled from HORIZONS."""

import os
from numpy import max
from skyfield import api
from skyfield.api import Topos
from skyfield.constants import AU_M

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

one_second = 1.0 / 24.0 / 60.0 / 60.0
arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
meter = 1.0 / AU_M

def ts():
    yield api.load.timescale()

def compare(value, expected_value, epsilon):
    if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
        assert max(abs(value - expected_value)) <= epsilon
    else:
        assert abs(value - expected_value) <= epsilon

def test_ecliptic_frame(ts):
    e = api.load('de421.bsp')
    jup = e['jupiter barycenter']
    astrometric = e['sun'].at(ts.utc(1980, 1, 1, 0, 0)).observe(jup)
    hlat, hlon, d = astrometric.ecliptic_latlon()
    compare(hlat.degrees, 1.013, 0.001)
    compare(hlon.degrees, 151.3229, 0.001)

def test_ecliptic_for_epoch_of_date(ts):
    e = api.load('de421.bsp')
    mars = e['mars barycenter']
    astrometric = e['earth'].at(ts.utc(1956, 1, 14, 6, 0, 0)).observe(mars)
    apparent = astrometric.apparent()
    hlat, hlon, d = apparent.ecliptic_latlon(epoch='date')
    compare(hlat.degrees, 0.4753402, 0.00001)
    compare(hlon.degrees, 240.0965633, 0.0002)

def test_ecliptic_for_epoch_of_date_array(ts):
    e = api.load('de421.bsp')
    sun = e['sun']
    astrometric = e['earth'].at(ts.utc(2005, 10, 1, [7, 8], 0, 0)).observe(sun)
    apparent = astrometric.apparent()
    hlat, hlon, d = apparent.ecliptic_latlon(epoch='date')
    compare(hlat.degrees[0], 0.0000488, 0.00001)
    compare(hlat.degrees[1], 0.0000474, 0.00001)
    compare(hlon.degrees[0], 188.2011122, 0.0002)
    compare(hlon.degrees[1], 188.2420983, 0.0002)

def test_fk4_frame(ts):
    e = api.load('de421.bsp')
    astrometric = e['earth'].at(ts.utc(1980, 1, 1, 0, 0)).observe(e['moon'])
    ra, dec, d = astrometric._to_spice_frame('B1950')
    print(ra._degrees, dec.degrees)
    compare(ra._degrees, 82.36186, 0.00006) # TODO: why is this not 0.00001?
    compare(dec.degrees, 18.53432, 0.00006)

def test_galactic_frame(ts):
    e = api.load('de421.bsp')
    astrometric = e['earth'].at(ts.utc(1980, 1, 1, 0, 0)).observe(e['moon'])
    glat, glon, d = astrometric.galactic_latlon()
    print(glat, glat.degrees, glon, glon.degrees)
    compare(glat.degrees, -8.047315, 0.005)  # TODO: awful! Track this down.
    compare(glon.degrees, 187.221794, 0.005)

def test_callisto_geometry(ts):
    e = api.load_file(_data_path('jup310-2053-10-08.bsp'))
    a = (e['callisto'] - e['earth']).at(ts.tdb(jd=2471184.5))
    print(a)
    compare(a.position.au,
      [-4.884815926454119E+00, -3.705745549073268E+00, -1.493487818022234E+00],
      0.001 * meter)
    compare(a.velocity.au_per_d,
      [9.604665478763035E-03, -1.552997751083403E-02, -6.678445860769302E-03],
      0.00001 * meter)

def test_callisto_astrometric(ts):
    e = api.load_file(_data_path('jup310-2053-10-08.bsp'))
    # This date was utc(2053, 10, 9), but new leap seconds keep breaking
    # the test, so:
    a = e['earth'].at(ts.tt(jd=2471184.5007775929)).observe(e['callisto'])
    ra, dec, distance = a.radec()
    compare(ra._degrees, 217.1839292, 0.001 * arcsecond)
    compare(dec.degrees, -13.6892791, 0.001 * arcsecond)
    compare(distance.au, 6.31079291776184, 0.1 * meter)

def test_boston_geometry():
    e = api.load_file(_data_path('jup310-2015-03-02.bsp'))
    t = api.load.timescale(delta_t=67.185390 + 0.5285957).tdb(2015, 3, 2)
    boston = e['earth'] + Topos((42, 21, 24.1), (-71, 3, 24.8),
                                x=0.003483, y=0.358609)
    a = (e['earth'] - boston).at(t)
    compare(a.position.km,
      [-1.764697476371664E+02, -4.717131288041386E+03, -4.274926422016179E+03],
      0.0027)  # TODO: try to get this < 1 meter

def test_moon_from_boston_geometry():
    e = api.load_file(_data_path('de430-2015-03-02.bsp'))
    t = api.load.timescale(delta_t=67.185390 + 0.5285957).tdb(2015, 3, 2)
    boston = e['earth'] + Topos((42, 21, 24.1), (-71, 3, 24.8),
                                x=0.003483, y=0.358609)
    a = (e['moon'] - boston).at(t)
    compare(a.position.au,
      [-1.341501206552443E-03, 2.190483327459023E-03, 6.839177007993498E-04],
      1.7 * meter)  # TODO: improve this

def test_moon_from_boston_astrometric():
    e = api.load_file(_data_path('de430-2015-03-02.bsp'))
    t = api.load.timescale(delta_t=67.185390 + 0.5285957).tdb(2015, 3, 2)
    boston = e['earth'] + Topos((42, 21, 24.1), (-71, 3, 24.8),
                                x=0.003483, y=0.358609)
    a = boston.at(t).observe(e['moon'])
    ra, dec, distance = a.radec()
    compare(ra._degrees, 121.4796470, 0.001 * arcsecond)
    compare(dec.degrees, 14.9108450, 0.001 * arcsecond)
    compare(distance.au, 0.00265828588792, 1.4 * meter)  # TODO: improve this
