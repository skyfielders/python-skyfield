"""Basic tests of the Skyfield API module and its contents."""

from assay import assert_raises
from skyfield import api
from skyfield import positionlib

def ts():
    yield api.load.timescale()

def test_whether_planets_have_radii():
    return # TODO: how will we support this again?
    assert api.mercury.radius.km == 2440.0
    for planet in api.nine_planets:
        assert planet.radius.km > 0.0

def test_sending_jd_that_is_not_a_julian_date():
    return # TODO: turn this back on, using one of the new ts method calls
    earth = api.load('de421.bsp')['earth']
    with assert_raises(ValueError, 'your "jd" argument is not a JulianDate: '):
        earth.at('blah')

def test_apparent_position_class(ts):
    e = api.load('de421.bsp')
    p = e['earth'].at(ts.utc((2014, 2, 9, 14, 50))).observe(e['mars']).apparent()
    assert isinstance(p, positionlib.Apparent)

def test_astrometric_position_class(ts):
    e = api.load('de421.bsp')
    p = e['earth'].at(ts.utc((2014, 2, 9, 14, 50))).observe(e['mars'])
    assert isinstance(p, positionlib.Astrometric)

def test_planet_position_class(ts):
    e = api.load('de421.bsp')
    p = e['mars'].at(ts.utc((2014, 2, 9, 14, 50)))
    assert isinstance(p, positionlib.Barycentric)

def test_star_position_class(ts):
    e = api.load('de421.bsp')
    star = api.Star(ra_hours=0, dec_degrees=0)
    p = e['earth'].at(ts.utc((2014, 2, 9, 15, 1))).observe(star)
    assert isinstance(p, positionlib.Astrometric)

def test_from_altaz_needs_topos():
    p = positionlib.ICRS([0.0, 0.0, 0.0])
    with assert_raises(ValueError, 'the orientation of the horizon'):
        p.from_altaz(alt_degrees=0, az_degrees=0)

def test_from_altaz_parameters(ts):
    e = api.load('de421.bsp')
    usno = e['earth'].topos('38.9215 N', '77.0669 W', elevation_m=92.0)
    jd = ts.tt(api.T0)
    p = usno.at(jd)
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
