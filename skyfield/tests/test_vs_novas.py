"""Compare the output of Skyfield with the same routines from NOVAS."""

import pytest
from numpy import array, einsum

from skyfield import (positionlib, earthlib, framelib, nutationlib,
                      jpllib, precessionlib, starlib, timescales)

from ..constants import ASEC2RAD, AU, AU_KM, DEG2RAD, T0
from ..functions import length_of
from ..timescales import JulianDate

# Since some users might run these tests without having installed our
# test dependencies, we detect import errors and skip these tests if the
# resources they need are not available.

try:
    import de405
    de405 = jpllib.Ephemeris(de405)
except ImportError:
    de405 = None

try:
    import novas
    import novas_de405
except ImportError:
    novas = None
else:
    import novas.compat as c
    import novas.compat.eph_manager

    jd_start, jd_end, number = c.eph_manager.ephem_open()  # needs novas_de405

    c_nutation = c.nutation
    import novas.compat.nutation  # overwrites nutation() function with module!

    TA = c.julian_date(1969, 7, 20, 20. + 18. / 60.)
    TB = c.julian_date(2012, 12, 21)
    TC = c.julian_date(2027, 8, 2, 10. + 7. / 60. + 50. / 3600.)

    D0 = 63.8285
    DA = 39.707
    DB = 66.8779
    DC = 72.  # http://maia.usno.navy.mil/ser7/deltat.preds

    P0 = (T0, D0)  # "pair 0"
    PA = (TA, DA)
    PB = (TB, DB)

arcminute = DEG2RAD / 60.0
arcsecond = arcminute / 60.0
arcsecond_in_hours = 24.0 / 360.0 / 60.0 / 60.0
arcsecond_in_degrees = 1.0 / 60.0 / 60.0
meter = 1.0 / AU

planet_codes = {
    'mercury': 1,
    'venus': 2,
    'mars': 4,
    'jupiter': 5,
    'saturn': 6,
    'uranus': 7,
    'neptune': 8,
    'pluto': 9,
    'sun': 10,
    'moon': 11,
    }

# Fixtures.

@pytest.fixture(params=[
    {'ut1': T0, 'delta_t': D0},
    {'ut1': TA, 'delta_t': DA},
    {'ut1': TB, 'delta_t': DB},
    {'ut1': TC, 'delta_t': DC},
    ])
def jd(request):
    # Build a new JulianDate each time, because some test cases need to
    # adjust the value of the date they are passed.
    return JulianDate(**request.param)

@pytest.fixture(params=[T0, TA, TB, TC])
def jd_float_or_vector(request):
    return request.param

@pytest.fixture(params=planet_codes.items())
def planet_name_and_code(request):
    return request.param

# Tests.

def eq(first, second, epsilon=None):
    """Test whether two floats are within `epsilon` of one another."""
    #print 'Significance of epsilon:', epsilon / second
    difference = abs(first - second)
    #print 'Difference relative to epsilon:', difference / epsilon
    if hasattr(first, 'shape') or hasattr(second, 'shape'):
        failed = difference.max() > epsilon
    else:
        failed = difference > epsilon
    if failed:
        appendix = ('\nbecause the difference is\n%r\ntimes too big'
                    % (abs(first - second) / epsilon)) if epsilon else ''
        raise AssertionError(
            '\n%r does not equal\n%r\nwithin the error bound\n%r%s'
            % (first, second, epsilon, appendix))


def test_star_deflected_by_jupiter(jd):
    star = c.make_cat_entry(
        star_name='Star', catalog='cat', star_num=101,
        ra=1.59132070233, dec=8.5958876464,
        pm_ra=0.0, pm_dec=0.0,
        parallax=0.0, rad_vel=0.0,
        )
    ra0, dec0 = c.app_star(jd.tt, star)

    earth = de405.earth
    star = starlib.Star(
        ra=1.59132070233, dec=8.5958876464,
        pm_ra=0.0, pm_dec=0.0,
        parallax=0.0, radial_velocity=0.0,
        )
    ra, dec, distance = earth(jd).observe(star).apparent().radec(epoch=jd)

    eq(ra0, ra.hours(), 1e-9 * arcsecond_in_hours)
    eq(dec0, dec.degrees(), 1e-9 * arcsecond_in_degrees)

# Tests of generating a full position or coordinate.

def test_astro_planet(jd, planet_name_and_code):
    planet_name, planet_code = planet_name_and_code

    obj = c.make_object(0, planet_code, 'planet', None)
    ra0, dec0, distance0 = c.astro_planet(jd.tt, obj)

    earth = de405.earth
    planet = getattr(de405, planet_name)
    e = earth(jd)
    distance = length_of((e - planet(jd)).position)
    ra, dec, d = e.observe(planet).radec()

    eq(ra0, ra.hours(), 1e-3 * arcsecond_in_hours)
    eq(dec0, dec.degrees(), 1e-3 * arcsecond_in_degrees)
    eq(distance0, distance, 0.5 * meter)

def test_virtual_planet(jd, planet_name_and_code):
    planet_name, planet_code = planet_name_and_code

    obj = c.make_object(0, planet_code, 'planet', None)
    ra0, dec0, distance0 = c.virtual_planet(jd.tt, obj)

    earth = de405.earth
    planet = getattr(de405, planet_name)
    e = earth(jd)
    distance = length_of((e - planet(jd)).position)
    ra, dec, d = e.observe(planet).apparent().radec()

    eq(ra0, ra.hours(), 0.001 * arcsecond_in_hours)
    eq(dec0, dec.degrees(), 0.001 * arcsecond_in_degrees)
    eq(distance0, distance, 0.5 * meter)

def test_app_planet(jd, planet_name_and_code):
    planet_name, planet_code = planet_name_and_code

    obj = c.make_object(0, planet_code, 'planet', None)
    ra0, dec0, distance0 = c.app_planet(jd.tt, obj)

    earth = de405.earth
    planet = getattr(de405, planet_name)
    e = earth(jd)
    distance = length_of((e - planet(jd)).position)
    ra, dec, d = e.observe(planet).apparent().radec(epoch=jd)

    eq(ra0, ra.hours(), 0.001 * arcsecond_in_hours)
    eq(dec0, dec.degrees(), 0.001 * arcsecond_in_degrees)
    eq(distance0, distance, 0.5 * meter)

def test_local_planet(jd, planet_name_and_code):
    position = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    ggr = positionlib.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.ephemeris = de405

    planet_name, planet_code = planet_name_and_code

    obj = c.make_object(0, planet_code, 'planet', None)
    ra0, dec0, distance0 = c.local_planet(jd.tt, jd.delta_t, obj, position)

    planet = getattr(de405, planet_name)
    g = ggr(jd)
    distance = length_of((g - planet(jd)).position)
    ra, dec, d = g.observe(planet).apparent().radec()

    eq(ra0, ra.hours(), 0.001 * arcsecond_in_hours)
    eq(dec0, dec.degrees(), 0.001 * arcsecond_in_degrees)
    eq(distance0, distance, 0.5 * meter)

def test_topo_planet(jd, planet_name_and_code):
    position = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    ggr = positionlib.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.ephemeris = de405

    planet_name, planet_code = planet_name_and_code

    obj = c.make_object(0, planet_code, 'planet', None)
    ra0, dec0, distance0 = c.topo_planet(jd.tt, jd.delta_t, obj, position)

    planet = getattr(de405, planet_name)
    g = ggr(jd)
    distance = length_of((g - planet(jd)).position)
    ra, dec, d = g.observe(planet).apparent().radec(epoch=jd)

    eq(ra0, ra.hours(), 0.001 * arcsecond_in_hours)
    eq(dec0, dec.degrees(), 0.001 * arcsecond_in_degrees)
    eq(distance0, distance, 0.5 * meter)

def test_altaz(jd, planet_name_and_code):
    """ Tests of generating a full position in altaz coordinates. Uses
        fixtures to iterate through date pairs and planets to generate
        individual tests.
    """
    planet_name, planet_code = planet_name_and_code
    position = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    ggr = positionlib.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.ephemeris = de405
    xp = yp = 0.0

    obj = c.make_object(0, planet_code, 'planet', None)
    ra, dec, dis = c.topo_planet(jd.tt, jd.delta_t, obj, position)
    (zd0, az0), (ra, dec) = c.equ2hor(jd.ut1, jd.delta_t, xp, yp,
                                      position, ra, dec)
    alt0 = 90.0 - zd0

    planet = getattr(de405, planet_name)
    g = ggr(jd)
    distance = length_of((g - planet(jd)).position)
    alt, az, d = g.observe(planet).apparent().altaz()

    eq(az0, az.degrees(), 0.001 * arcsecond_in_degrees)
    eq(alt0, alt.degrees(), 0.001 * arcsecond_in_degrees)
    eq(dis, distance, 0.5 * meter)

# Tests for Basic Functions

def test_cal_date():
    for jd in 0.0, 2414988.5, 2415020.31352, 2442249.5, 2456335.2428472:
        assert c.cal_date(jd) == timescales.cal_date(jd)

def test_earth_rotation_angle(jd_float_or_vector):
    jd_ut1 = jd_float_or_vector
    u = c.era(jd_ut1)
    v = earthlib.earth_rotation_angle(jd_ut1)
    epsilon = 1e-12  # degrees; 14 to 15 digits of agreement
    eq(u, v, epsilon)

def test_earth_tilt(jd):
    u = c.e_tilt(jd.tdb)
    v = nutationlib.earth_tilt(jd)
    epsilon = 1e-9  # 9 to 11 digits of agreement; why not more?
    eq(array(u), array(v), epsilon)

def test_equation_of_the_equinoxes_complimentary_terms(jd_float_or_vector):
    jd_tt = jd_float_or_vector
    u = c.ee_ct(jd_tt, 0.0, 0)
    v = nutationlib.equation_of_the_equinoxes_complimentary_terms(jd_tt)

    epsilon = 1e-22  # radians; 14 digits of agreement
    eq(u, v, epsilon)

def test_frame_tie():
    xyz = array([1.1, 1.2, 1.3])
    epsilon = 0.0  # perfect
    eq(c.frame_tie(xyz, 0), xyz.dot(framelib.ICRS_to_J2000), epsilon)
    eq(c.frame_tie(xyz, -1), xyz.dot(framelib.J2000_to_ICRS), epsilon)

def test_fundamental_arguments(jd_float_or_vector):
    jd_tdb = jd_float_or_vector
    t = jcentury(jd_tdb)
    u = c.fund_args(t)
    v = nutationlib.fundamental_arguments(t)

    epsilon = 1e-12  # radians; 13 digits of agreement
    eq(u, v, epsilon)

def test_geocentric_position_and_velocity(jd):
    observer = c.make_observer_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    posu, velu = c.geo_posvel(jd.tt, jd.delta_t, observer)

    topos = positionlib.Topos('75 W', '45 N', elevation=0.0,
                              temperature=10.0, pressure=1010.0)
    posv, velv = earthlib.geocentric_position_and_velocity(topos, jd)

    epsilon = 1e-6 * meter  # 13 to 14 digits of agreement

    eq(posu, posv, epsilon)
    eq(velu, velv, epsilon)

def test_iau2000a(jd_float_or_vector):
    jd_tt = jd_float_or_vector
    psi0, eps0 = c.nutation.iau2000a(jd_tt, 0.0)
    psi1, eps1 = nutationlib.iau2000a(jd_tt)
    to_tenths_of_microarcseconds = 1e7 / ASEC2RAD

    epsilon = 4e-6  # tenths of micro arcseconds; 13 digits of precision

    eq(psi0 * to_tenths_of_microarcseconds, psi1, epsilon)
    eq(eps0 * to_tenths_of_microarcseconds, eps1, epsilon)

def test_julian_date():
    epsilon = 0.0  # perfect
    for args in (
          (-4712, 1, 1, 0.0),
          (-4712, 3, 1, 0.0),
          (-4712, 12, 31, 0.5),
          (-241, 3, 25, 19.0),
          (530, 9, 27, 23.5),
          (1976, 3, 7, 12.5),
          (2000, 1, 1, 0.0),
          ):
        eq(c.julian_date(*args), timescales.julian_date(*args), epsilon)

def test_mean_obliq(jd_float_or_vector):
    jd_tdb = jd_float_or_vector
    u = c.mean_obliq(jd_tdb)
    v = nutationlib.mean_obliquity(jd_tdb)
    epsilon = 0.0  # perfect
    eq(u, v, epsilon)

def test_nutation(jd):
    xyz = [1.1, 1.2, 1.3]
    u = c_nutation(jd.tdb, xyz)
    xyz = array(xyz)
    v = einsum('i...,ij...->j...', xyz, nutationlib.compute_nutation(jd))
    epsilon = 1e-14  # 14 digits of agreement
    eq(u, v, epsilon)

def test_precession(jd_float_or_vector):
    jd_tdb = jd_float_or_vector
    xyz = [1.1, 1.2, 1.3]
    u = c.precession(T0, xyz, jd_tdb)
    matrix_or_matrices = precessionlib.compute_precession(jd_tdb)
    v = einsum('ij...,i...->j...', matrix_or_matrices, array(xyz))
    epsilon = 1e-15  # 15 digits of agreement
    eq(u, v, epsilon)

def test_sidereal_time_with_zero_delta_t(jd):
    jd.delta_t = 0.0
    u = c.sidereal_time(jd.ut1, 0.0, 0.0, False, True)
    v = earthlib.sidereal_time(jd)
    epsilon = 1e-13  # days; 14 digits of agreement
    eq(u, v, epsilon)

def test_sidereal_time_with_nonzero_delta_t(jd):
    u = c.sidereal_time(jd.ut1, 0.0, jd.delta_t, False, True)
    v = earthlib.sidereal_time(jd)
    epsilon = 1e-13  # days; 14 digits of agreement
    eq(u, v, epsilon)

def test_starvectors():
    p, v = c.starvectors(c.make_cat_entry(
            'POLARIS', 'HIP', 0, 2.530301028, 89.264109444,
            44.22, -11.75, 7.56, -17.4))

    star = starlib.Star(2.530301028, 89.264109444,
                        44.22, -11.75, 7.56, -17.4)

    p_epsilon = 1e-10  # AU; 16 digits of agreement
    v_epsilon = 1e-17  # AU/day; 15 digits of agreement

    eq(p, star._position, p_epsilon)
    eq(v, star._velocity, v_epsilon)

def test_ter2cel(jd):
    jd_low = 0.0
    xp = yp = 0.0

    position = array([1.1, 1.2, 1.3])

    theirs = c.ter2cel(jd.ut1, jd_low, jd.delta_t, xp, yp, position)
    ours = positionlib.ITRF_to_GCRS(jd, position)

    epsilon = 1e-13  # 13 digits of agreement
    eq(theirs, ours, epsilon)

def test_terra():
    observer = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)

    # Note that this class stands in for a NOVAS Topos structure, but
    # not for our own Topos class!
    class Topos(object):
        latitude = 45.0 * DEG2RAD
        longitude = -75.0 * DEG2RAD
        elevation = 0.0
    topos = Topos()

    pos0, vel0 = array(c.terra(observer, 11.0))
    pos1, vel1 = array(c.terra(observer, 23.9))

    posn, veln = earthlib.terra(topos, array([11.0, 23.9]))

    epsilon = 1e-8 * meter  # 14 digits of agreement

    eq(pos0, posn[:,0], epsilon)
    eq(pos1, posn[:,1], epsilon)
    eq(vel0, veln[:,0], epsilon)
    eq(vel1, veln[:,1], epsilon)

def test_tdb2tt(jd_float_or_vector):
    jd_tdb = jd_float_or_vector
    u = c.tdb2tt(jd_tdb)[1]
    v = timescales.tdb_minus_tt(jd_tdb)
    epsilon_seconds = 1e-16  # 11 or 12 digits of agreement; why not more?
    eq(u, v, epsilon_seconds)

def jcentury(t):
    return (t - T0) / 36525.0
