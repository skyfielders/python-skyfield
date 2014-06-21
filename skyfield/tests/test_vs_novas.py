"""Compare the output of Skyfield with the same routines from NOVAS."""

import pytest
from numpy import array, einsum

from skyfield import (positionlib, earthlib, framelib, nutationlib,
                      jpllib, precessionlib, starlib, timelib)

from ..constants import ASEC2RAD, AU_M, DEG2RAD, T0
from ..functions import length_of
from ..timelib import JulianDate

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
meter = 1.0 / AU_M

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
    {'tt': T0, 'delta_t': D0},
    {'tt': TA, 'delta_t': DA},
    {'tt': TB, 'delta_t': DB},
    {'tt': TC, 'delta_t': DC},
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
    # if epsilon:
    #     print 'Difference relative to epsilon:', difference / epsilon
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
        ra_hours=1.59132070233, dec_degrees=8.5958876464,
        ra_mas_per_year=0.0, dec_mas_per_year=0.0,
        parallax_mas=0.0, radial_km_per_s=0.0,
        )
    ra, dec, distance = earth(jd).observe(star).apparent().radec(epoch=jd)

    eq(ra0, ra.hours, 1e-9 * arcsecond_in_hours)
    eq(dec0, dec.degrees, 1e-9 * arcsecond_in_degrees)

# Tests for Basic Functions

def test_geocentric_position_and_velocity(jd):
    observer = c.make_observer_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
    posu, velu = c.geo_posvel(jd.tt, jd.delta_t, observer)

    topos = positionlib.Topos(latitude_degrees=45.0, longitude_degrees=-75.0,
                              elevation_m=0.0)
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
        eq(c.julian_date(*args), timelib.julian_date(*args), epsilon)

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
    v = einsum('ij...,j...->i...', nutationlib.compute_nutation(jd), xyz)
    epsilon = 1e-14  # 14 digits of agreement
    eq(u, v, epsilon)

def test_precession(jd_float_or_vector):
    jd_tdb = jd_float_or_vector
    xyz = [1.1, 1.2, 1.3]
    u = c.precession(T0, xyz, jd_tdb)
    matrix_or_matrices = precessionlib.compute_precession(jd_tdb)
    v = einsum('ij...,j...->i...', matrix_or_matrices, array(xyz))
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

    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    p_epsilon = 1e3 * meter   # not bad for something 27 million AU distant
    v_epsilon = 1e-6 * meter  # per day

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
    topos = positionlib.Topos(latitude_degrees=45.0, longitude_degrees=-75.0,
                              elevation_m=0.0)

    pos0, vel0 = array(c.terra(observer, 11.0))
    pos1, vel1 = array(c.terra(observer, 23.9))

    posn, veln = earthlib.terra(topos, array([11.0, 23.9]))

    epsilon = 1e-7 * meter  # 13 digits of agreement

    eq(pos0, posn[:,0], epsilon)
    eq(pos1, posn[:,1], epsilon)
    eq(vel0, veln[:,0], epsilon)
    eq(vel1, veln[:,1], epsilon)

def test_tdb2tt(jd_float_or_vector):
    jd_tdb = jd_float_or_vector
    u = c.tdb2tt(jd_tdb)[1]
    v = timelib.tdb_minus_tt(jd_tdb)
    epsilon_seconds = 1e-16  # 11 or 12 digits of agreement; why not more?
    eq(u, v, epsilon_seconds)

def jcentury(t):
    return (t - T0) / 36525.0
