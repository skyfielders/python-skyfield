"""Compare the output of Skyfield with the same routines from NOVAS."""

import pytest
from numpy import array, einsum, rollaxis

from skyfield import (coordinates, earthlib, framelib, nutationlib,
                      jpllib, precessionlib, starlib, timescales)

from ..constants import T0, DEG2RAD, AU_KM, TAU
from ..timescales import JulianDate

# Since some users might run these tests without having installed our
# test dependencies, we detect import errors and skip these tests if the
# resources they need are not available.

try:
    import de405
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
meter = 1.0 / AU_KM

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

planets_to_test = planet_codes.keys()

# Fixtures.

class Box(object):
    """Hides values, like NumPy arrays, from py.test fixture bugs."""
    def __init__(self, value):
        self.boxed_value = value

jd_vector = array([T0, TA, TB, TC])
jd_vector.flags.writeable = False

@pytest.fixture(params=[
    Box({'ut1': T0, 'delta_t': D0}),
    Box({'ut1': TA, 'delta_t': DA}),
    Box({'ut1': TB, 'delta_t': DB}),
    Box({'ut1': TC, 'delta_t': DC}),
    Box({'ut1': [T0, TA, TB, TC], 'delta_t': [D0, DA, DB, DC]}),
    ])
def jd(request):
    # Build a new JulianDate each time, because some test cases need to
    # adjust the value of the date they are passed.
    return JulianDate(**request.param.boxed_value)

@pytest.fixture(params=[Box(T0), Box(TA), Box(TB), Box(TC), Box(jd_vector)])
def jd_float_or_vector(request):
    return request.param

@pytest.fixture(params=[P0, PA, PB])
def timepairs(request):
    return request.param

@pytest.fixture(params=planet_codes.items())
def planets_list(request):
    return request.param

# Helpers.

def vcall(function, *args):
    """Call a function once or many times, on whether any args are arrays.

    If no arguments are NumPy arrays, then `function` is called once
    with the arguments and its return value is returned.  Otherwise,
    `function` is called as many times as there are values in the
    arrays, and an array of return values is returned.

    """
    lengths = [arg.shape[0] for arg in args if getattr(arg, 'shape', ())]
    if not lengths:
        return function(*args)
    length = min(lengths)
    arglists = [list() for i in range(length)]
    for arg in args:
        if hasattr(arg, 'shape'):
            for i, arglist in enumerate(arglists):
                arglist.append(arg[i])
        else:
            for arglist in arglists:
                arglist.append(arg)
    result = array([function(*arglist) for arglist in arglists])
    return rollaxis(result, 0, len(result.shape))

# Tests.

emp = jpllib.Ephemeris(de405)

def eq(first, second, epsilon=None):
    """Test whether two floats are within `epsilon` of one another."""
    if hasattr(first, 'shape') or hasattr(second, 'shape'):
        failed = abs(first - second).max() > epsilon
    else:
        failed = abs(first - second) > epsilon
    if failed:
        appendix = ('\nbecause the difference is\n%r\ntimes too big'
                    % (abs(first - second) / epsilon)) if epsilon else ''
        raise AssertionError(
            '%r\ndoes not equal\n%r\nwithin the error bound\n%r%s'
            % (first, second, epsilon, appendix))


def test_new_star_deflected_by_jupiter(timepairs):
    """ Tests of generating a stellar position. """
    jd_tt = timepairs[0]
    star = c.make_cat_entry(
        star_name='Star', catalog='cat', star_num=101,
        ra=1.59132070233, dec=8.5958876464,
        pm_ra=0.0, pm_dec=0.0,
        parallax=0.0, rad_vel=0.0,
        )
    ra, dec = c.app_star(jd_tt, star)

    earth = emp.earth
    star = starlib.Star(
        ra=1.59132070233, dec=8.5958876464,
        pm_ra=0.0, pm_dec=0.0,
        parallax=0.0, radial_velocity=0.0,
        )
    jd = JulianDate(tt=jd_tt)
    g = star.observe_from(earth(jd)).apparent()

    eq(ra * TAU / 24.0, g.ra, 0.001 * arcsecond)
    eq(dec * TAU / 360.0, g.dec, 0.001 * arcsecond)

# Tests of generating a full position or coordinate.

def test_astro_planet(timepairs, planets_list):
    jd_tt = timepairs[0]
    planet_name = planets_list[0]
    planet_code = planets_list[1]

    obj = c.make_object(0, planet_code, 'planet', None)
    ra, dec, dis = c.astro_planet(jd_tt, obj)

    earth = emp.earth
    planet = getattr(emp, planet_name)
    jd = JulianDate(tt=jd_tt)
    g = planet.observe_from(earth(jd)).astrometric()

    eq(ra * TAU / 24.0, g.ra, 0.001 * arcsecond)
    eq(dec * TAU / 360.0, g.dec, 0.001 * arcsecond)
    eq(dis, g.distance, 0.001 * meter)


def test_app_planet(timepairs, planets_list):
    jd_tt = timepairs[0]
    planet_name = planets_list[0]
    planet_code = planets_list[1]

    obj = c.make_object(0, planet_code, 'planet', None)
    ra, dec, dis = c.app_planet(jd_tt, obj)

    earth = emp.earth
    planet = getattr(emp, planet_name)
    jd = JulianDate(tt=jd_tt)
    g = planet.observe_from(earth(jd)).apparent()

    eq(ra * TAU / 24.0, g.ra, 0.001 * arcsecond)
    eq(dec * TAU / 360.0, g.dec, 0.001 * arcsecond)
    eq(dis, g.distance, 0.001 * meter)

def test_topo_planet(timepairs, planets_list):
    position = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    ggr = coordinates.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.ephemeris = emp

    jd_tt = timepairs[0]
    delta_t = timepairs[1]
    planet_name = planets_list[0]
    planet_code = planets_list[1]

    obj = c.make_object(0, planet_code, 'planet', None)
    ra, dec, dis = c.topo_planet(jd_tt, delta_t, obj, position)

    planet = getattr(emp, planet_name)
    jd = JulianDate(tt=jd_tt, delta_t=delta_t)
    g = ggr(jd).observe(planet).apparent()

    eq(ra * TAU / 24.0, g.ra, 0.001 * arcsecond)
    eq(dec * TAU / 360.0, g.dec, 0.001 * arcsecond)
    eq(dis, g.distance, 0.001 * meter)


def test_new_horizontal(timepairs, planets_list):
    """ Tests of generating a full position in horizontal coordinates. Uses
        fixtures to iterate through date pairs and planets to generate
        individual tests.
    """
    jd_tt = timepairs[0]
    delta_t = timepairs[1]
    planet_name = planets_list[0]
    planet_code = planets_list[1]
    position = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    ggr = coordinates.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.ephemeris = emp
    xp = yp = 0.0

    # replaces the for loop
    obj = c.make_object(0, planet_code, 'planet', None)
    ra, dec, dis = c.topo_planet(jd_tt, delta_t, obj, position)
    jd_ut1 = jd_tt - delta_t / 86400.0
    (zd, az), (ra, dec) = c.equ2hor(
        jd_ut1, delta_t, xp, yp, position, ra, dec, ref_option=0)
    planet = getattr(emp, planet_name)
    jd = JulianDate(tt=jd_tt, delta_t=delta_t)
    h = ggr(jd).observe(planet).apparent().horizontal()

    eq(zd * TAU / 360.0, h.zd, 0.001 * arcsecond)
    eq(az * TAU / 360.0, h.az, 0.001 * arcsecond)
    eq(0.25 * TAU - zd * TAU / 360.0, h.alt, 0.001 * arcsecond)
    eq(dis, h.distance, 0.001 * meter)

# Tests for Basic Functions

def test_cal_date():
    for jd in 0.0, 2414988.5, 2415020.31352, 2442249.5, 2456335.2428472:
        assert c.cal_date(jd) == timescales.cal_date(jd)

def test_earth_rotation_angle(jd_float_or_vector):
    epsilon = 1e-12
    jd_ut1 = jd_float_or_vector.boxed_value
    u = vcall(c.era, jd_ut1)
    v = earthlib.earth_rotation_angle(jd_ut1)
    eq(u, v, epsilon)

def test_earth_tilt(jd):
    epsilon = 1e-9
    u = vcall(c.e_tilt, jd.tdb)
    v = nutationlib.earth_tilt(jd)
    eq(u, array(v), epsilon)

def test_equation_of_the_equinoxes_complimentary_terms(jd_float_or_vector):
    epsilon = 1e-22
    jd_tt = jd_float_or_vector.boxed_value
    u = vcall(c.ee_ct, jd_tt, 0.0, 0)
    v = nutationlib.equation_of_the_equinoxes_complimentary_terms(jd_tt)
    eq(u, v, epsilon)

def test_frame_tie():
    epsilon = 1e-15
    xyz = array([1.0, 2.0, 3.0])

    eq(c.frame_tie(xyz, 0), xyz.dot(framelib.ICRS_to_J2000), epsilon)
    eq(c.frame_tie(xyz, -1), xyz.dot(framelib.J2000_to_ICRS), epsilon)

def test_fundamental_arguments(jd_float_or_vector):
    epsilon = 1e-12
    jd_tdb = jd_float_or_vector.boxed_value
    t = jcentury(jd_tdb)
    u = vcall(c.fund_args, t)
    v = nutationlib.fundamental_arguments(t)
    eq(u, v, epsilon)

def test_geocentric_position_and_velocity(jd):
    epsilon = 1e-13
    jd.delta_t = delta_t = 0.0  # TODO: relax this limitation?

    observer = c.make_observer_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    posu, velu = vcall(c.geo_posvel, jd.tt, delta_t, observer)

    topos = coordinates.Topos('75 W', '45 N', elevation=0.0,
                              temperature=10.0, pressure=1010.0)
    posv, velv = earthlib.geocentric_position_and_velocity(topos, jd)

    eq(posu, posv, epsilon)
    eq(velu, velv, epsilon)

def test_iau2000a(jd_float_or_vector):
    epsilon = 2e-18
    jd_tt = jd_float_or_vector.boxed_value
    psi0, eps0 = vcall(c.nutation.iau2000a, jd_tt, 0.0)
    psi1, eps1 = nutationlib.iau2000a(jd_tt)
    eq(psi0, psi1, epsilon)
    eq(eps0, eps1, epsilon)

def test_julian_date():
    epsilon = 0.0
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
    epsilon = 0
    jd_tdb = jd_float_or_vector.boxed_value
    u = vcall(c.mean_obliq, jd_tdb)
    v = nutationlib.mean_obliquity(jd_tdb)
    eq(u, v, epsilon)

def test_nutation(jd):
    epsilon = 1e-15
    xyz = [1.0, 2.0, 3.0]
    u = vcall(c_nutation, jd.tt, xyz)  # TODO: shouldn't this be jd.tdb?
    xyz = array(xyz)
    v = einsum('i...,ij...->j...', xyz, nutationlib.compute_nutation(jd))
    eq(u, v, epsilon)

def test_precession(jd_float_or_vector):
    epsilon = 1e-15
    jd_tdb = jd_float_or_vector.boxed_value
    xyz = [1.0, 2.0, 3.0]
    u = vcall(c.precession, T0, xyz, jd_tdb)
    xyz = array(xyz)
    matrix_or_matrices = precessionlib.compute_precession(jd_tdb)
    v = einsum('i...,ij...->j...', array(xyz), matrix_or_matrices)
    eq(u, v, epsilon)

def test_sidereal_time_with_zero_delta_t(jd):
    epsilon = 1e-13
    jd.delta_t = 0.0
    u = vcall(c.sidereal_time, jd.ut1, 0.0, 0.0, False, True)
    v = earthlib.sidereal_time(jd)
    eq(u, v, epsilon)

def test_sidereal_time_with_nonzero_delta_t(jd):
    epsilon = 1e-13
    u = vcall(c.sidereal_time, jd.ut1, 0.0, jd.delta_t, False, True)
    v = earthlib.sidereal_time(jd)
    eq(u, v, epsilon)

def test_starvectors():
    epsilon = 1e-10

    p, v = c.starvectors(c.make_cat_entry(
            'POLARIS', 'HIP', 0, 2.530301028, 89.264109444,
            44.22, -11.75, 7.56, -17.4))

    star = starlib.Star(2.530301028, 89.264109444,
                        44.22, -11.75, 7.56, -17.4)

    eq(p, star._position.reshape(3), epsilon)
    eq(v, star._velocity.reshape(3), epsilon)

def test_terra():
    epsilon = 1e-18

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

    eq(pos0, posn[:,0], epsilon)
    eq(pos1, posn[:,1], epsilon)
    eq(vel0, veln[:,0], epsilon)
    eq(vel1, veln[:,1], epsilon)

def test_tdb2tt(jd_float_or_vector):
    epsilon = 1e-16
    jd_tdb = jd_float_or_vector.boxed_value
    u = vcall(c.tdb2tt, jd_tdb)[1]
    v = timescales.tdb_minus_tt(jd_tdb)
    eq(u, v, epsilon)

def jcentury(t):
    return (t - T0) / 36525.0
