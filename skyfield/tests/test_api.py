"""Tests of the Skyfield `api` module and user-facing exceptions."""

import numpy as np
from assay import assert_raises
from skyfield import api, positionlib
from skyfield.api import Topos
from skyfield.errors import EphemerisRangeError

def ts():
    yield api.load.timescale()

def test_sending_jd_that_is_not_a_julian_date():
    earth = api.load('de421.bsp')['earth']
    with assert_raises(ValueError, r"please provide the at\(\) method"
                       " with a Time instance as its argument,"
                       " instead of the value 'blah'"):
        earth.at('blah')

def test_observe_and_apparent_survive_an_empty_time_array(ts):
    e = api.load('de421.bsp')
    e['earth'].at(ts.utc(2024, 2, 4, 7, [])).observe(e['mars']).apparent()

def test_apparent_position_class(ts):
    e = api.load('de421.bsp')
    p = e['earth'].at(ts.utc(2014, 2, 9, 14, 50)).observe(e['mars']).apparent()
    assert isinstance(p, positionlib.Apparent)

def test_astrometric_position_class(ts):
    e = api.load('de421.bsp')
    p = e['earth'].at(ts.utc(2014, 2, 9, 14, 50)).observe(e['mars'])
    assert isinstance(p, positionlib.Astrometric)

def test_astrometric_position_does_not_allow_altaz(ts):
    e = api.load('de421.bsp')
    o = e['earth'] + api.wgs84.latlon(36.7138, -112.2169)
    a = o.at(ts.utc(2014, 2, 9, 14, 50)).observe(e['mars'])
    with assert_raises(ValueError, 'it is not useful to call .altaz'):
        a.altaz()

def test_ephemeris_contains_method(ts):
    e = api.load('de421.bsp')
    assert (399 in e) is True
    assert (398 in e) is False
    assert ('earth' in e) is True
    assert ('Earth' in e) is True
    assert ('EARTH' in e) is True
    assert ('ceres' in e) is False

def test_exception_raised_for_dates_outside_ephemeris(ts):
    eph = api.load('de421.bsp')
    message = (
        'ephemeris segment only covers dates 1899-07-29 through 2053-10-09'
    )
    with assert_raises(EphemerisRangeError, message) as a:
        eph['earth'].at(ts.tt(4096))

    e = a.exception
    assert e.args == (message,)
    assert e.start_time.tdb == 2414864.5
    assert e.end_time.tdb == 2471184.5
    assert e.time_mask == [True]
    assert e.segment is eph['earth'].vector_functions[0].spk_segment

def test_planet_position_class(ts):
    e = api.load('de421.bsp')
    p = e['mars'].at(ts.utc(2014, 2, 9, 14, 50))
    assert isinstance(p, positionlib.Barycentric)

def test_star_position_class(ts):
    e = api.load('de421.bsp')
    star = api.Star(ra_hours=0, dec_degrees=0)
    p = e['earth'].at(ts.utc(2014, 2, 9, 15, 1)).observe(star)
    assert isinstance(p, positionlib.Astrometric)

def test_star_vector_from_earth(ts):
    t = ts.tt_jd(api.T0)
    eph = api.load('de421.bsp')
    e = eph['earth'].at(t)

    star = api.Star(ra_hours=[1.0, 2.0], dec_degrees=[+3.0, +4.0])
    p = e.observe(star)
    assert p.xyz.au.shape == (3, 2)
    assert p.velocity.au_per_d.shape == (3, 2)
    assert p.t.shape == (2,)
    assert (p.t.tt == api.T0).all()
    a = p.apparent()

    a1 = e.observe(api.Star(ra_hours=1.0, dec_degrees=+3.0)).apparent()
    a2 = e.observe(api.Star(ra_hours=2.0, dec_degrees=+4.0)).apparent()
    assert (a1.xyz.au == a.xyz.au[:,0]).all()
    assert (a2.xyz.au == a.xyz.au[:,1]).all()

def test_star_vector_from_topos(ts):
    t = ts.tt_jd(api.T0)
    eph = api.load('de421.bsp')
    boston = eph['earth'] + Topos('42.3583 N', '71.0636 W')
    b = boston.at(t)

    star = api.Star(ra_hours=[1.0, 2.0], dec_degrees=[+3.0, +4.0])
    p = b.observe(star)
    assert p.xyz.au.shape == (3, 2)
    assert p.velocity.au_per_d.shape == (3, 2)
    assert p.t.shape == (2,)
    assert (p.t.tt == api.T0).all()
    a = p.apparent()

    a1 = b.observe(api.Star(ra_hours=1.0, dec_degrees=+3.0)).apparent()
    a2 = b.observe(api.Star(ra_hours=2.0, dec_degrees=+4.0)).apparent()
    assert (a1.xyz.au == a.xyz.au[:,0]).all()
    assert (a2.xyz.au == a.xyz.au[:,1]).all()

def test_hadec_needs_a_longitude(ts):
    e = api.load('de421.bsp')
    earth = e['earth']
    moon = e['moon']
    apparent = earth.at(ts.utc(2016)).observe(moon)
    with assert_raises(ValueError, 'from a specific latitude and longitude'):
        apparent.hadec()

def test_altaz_needs_topos(ts):
    e = api.load('de421.bsp')
    earth = e['earth']
    moon = e['moon']
    apparent = earth.at(ts.utc(2016)).observe(moon).apparent()
    with assert_raises(ValueError, 'from a specific latitude and longitude'):
        apparent.altaz()

def test_from_altaz_needs_topos():
    p = positionlib.ICRF([0.0, 0.0, 0.0])
    with assert_raises(ValueError, 'to compute altitude and azimuth'):
        p.from_altaz(alt_degrees=0, az_degrees=0)

def test_from_altaz_parameters(ts):
    usno = Topos('38.9215 N', '77.0669 W', elevation_m=92.0)
    t = ts.tt(jd=api.T0)
    p = usno.at(t)
    a = api.Angle(degrees=10.0)
    d = api.Distance(au=0.234)
    with assert_raises(ValueError, 'the alt= parameter with an Angle'):
        p.from_altaz(alt='Bad value', alt_degrees=0, az_degrees=0)
    with assert_raises(ValueError, 'the az= parameter with an Angle'):
        p.from_altaz(az='Bad value', alt_degrees=0, az_degrees=0)
    p.from_altaz(alt=a, alt_degrees='bad', az_degrees=0)
    p.from_altaz(az=a, alt_degrees=0, az_degrees='bad')
    assert str(p.from_altaz(alt=a, az=a).distance()) == '0.1 au'
    assert str(p.from_altaz(alt=a, az=a, distance=d).distance()) == '0.234 au'

def test_github_81(ts):
    t = ts.utc(1980, 1, 1)
    assert t == t
    assert (t == 61) is False  # used to die with AttributeError for "tt"

def test_github_500_does_zero_position_trigger_numpy_warnings(ts):
    zero_vector = np.array([0.0, 0.0, 0.0])
    t = ts.utc(2020, 12, 14)
    p = positionlib.Astrometric(zero_vector, zero_vector, t=t, center=0)
    p._ephemeris = api.load('de421.bsp')
    with np.errstate(all='raise'):
        p.apparent().radec()
