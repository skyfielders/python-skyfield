"""Basic tests of the Skyfield API module and its contents."""

from assay import assert_raises
from skyfield import api, positionlib
from skyfield.api import Topos

def ts():
    yield api.load.timescale()

def test_whether_planets_have_radii():
    return # TODO: how will we support this again?
    assert api.mercury.radius.km == 2440.0
    for planet in api.nine_planets:
        assert planet.radius.km > 0.0

def test_sending_jd_that_is_not_a_julian_date():
    earth = api.load('de421.bsp')['earth']
    with assert_raises(ValueError, r"please provide the at\(\) method"
                       " with a Time instance as its argument,"
                       " instead of the value 'blah'"):
        earth.at('blah')

def test_apparent_position_class(ts):
    e = api.load('de421.bsp')
    p = e['earth'].at(ts.utc(2014, 2, 9, 14, 50)).observe(e['mars']).apparent()
    assert isinstance(p, positionlib.Apparent)

def test_astrometric_position_class(ts):
    e = api.load('de421.bsp')
    p = e['earth'].at(ts.utc(2014, 2, 9, 14, 50)).observe(e['mars'])
    assert isinstance(p, positionlib.Astrometric)

def test_ephemeris_contains_method(ts):
    e = api.load('de421.bsp')
    assert (399 in e) is True
    assert (398 in e) is False
    assert ('earth' in e) is True
    assert ('Earth' in e) is True
    assert ('EARTH' in e) is True
    assert ('ceres' in e) is False

def test_planet_position_class(ts):
    e = api.load('de421.bsp')
    p = e['mars'].at(ts.utc(2014, 2, 9, 14, 50))
    assert isinstance(p, positionlib.Barycentric)

def test_star_position_class(ts):
    e = api.load('de421.bsp')
    star = api.Star(ra_hours=0, dec_degrees=0)
    p = e['earth'].at(ts.utc(2014, 2, 9, 15, 1)).observe(star)
    assert isinstance(p, positionlib.Astrometric)

def test_altaz_needs_topos(ts):
    e = api.load('de421.bsp')
    earth = e['earth']
    moon = e['moon']
    with assert_raises(ValueError, 'using a Topos instance'):
        earth.at(ts.utc(2016)).observe(moon).apparent().altaz()

def test_from_altaz_needs_topos():
    p = positionlib.ICRF([0.0, 0.0, 0.0])
    with assert_raises(ValueError, 'the orientation of the horizon'):
        p.from_altaz(alt_degrees=0, az_degrees=0)

def test_from_altaz_parameters(ts):
    e = api.load('de421.bsp')
    usno = e['earth'] + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)
    t = ts.tt(jd=api.T0)
    p = usno.at(t)
    a = api.Angle(degrees=10.0)
    with assert_raises(ValueError, 'the alt= parameter with an Angle'):
        p.from_altaz(alt='Bad value', alt_degrees=0, az_degrees=0)
    with assert_raises(ValueError, 'the az= parameter with an Angle'):
        p.from_altaz(az='Bad value', alt_degrees=0, az_degrees=0)
    p.from_altaz(alt=a, alt_degrees='bad', az_degrees=0)
    p.from_altaz(az=a, alt_degrees=0, az_degrees='bad')

def test_named_star_throws_valueerror():
    with assert_raises(ValueError, 'No star named foo known to skyfield'):
        api.NamedStar('foo')

def test_named_star_returns_star():
    s = api.NamedStar('Polaris')
    assert isinstance(s, api.Star)

def test_github_81(ts):
    t = ts.utc(1980, 1, 1)
    assert t == t
    assert (t == 61) is False  # used to die with AttributeError for "tt"
