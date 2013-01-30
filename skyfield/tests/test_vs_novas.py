"""Compare the output of Skyfield with the same routines from NOVAS."""

from itertools import product
from numpy import array
from unittest import TestCase

from skyfield import (angles, coordinates, earthlib, framelib, nutationlib,
                      planets, precessionlib, timescales)

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

    TA = c.julian_date(1969, 7, 20, 20. + 18./60.)  # arbitrary test date
    TB = c.julian_date(2012, 12, 21)                # arbitrary test date

tau = angles.tau
degree = tau / 360.0
arcminute = degree / 60.0
arcsecond = arcminute / 60.0
meter = 1.0 / earthlib.AU_KM
T0 = timescales.T0

A0 = array([T0])
AA = array([TA])
AB = array([TB])

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

class NOVASTests(TestCase):

    delta = 'the delta needs to be specified at the top of each test'

    @classmethod
    def setUpClass(cls):
        if de405 is None or novas is None:
            cls.__unittest_skip__ = True
            return
        cls.e = planets.Ephemeris(de405)

    def eq(self, first, second, delta=None):
        if delta is None:
            delta = self.delta
        if abs(first - second) > delta:
            raise AssertionError(
                '%r != %r within %r because the difference is %r times too big'
                % (first, second, delta, abs(first - second) / delta))

    # Tests of generating a full position or coordinate.

    def test_astro_planet(self):

        for t, name in product((T0, TA, TB), planets_to_test):
            obj = c.make_object(0, planet_codes[name], b'planet', None)
            ra, dec, dis = c.astro_planet(t, obj)

            earth = self.e.earth
            planet = getattr(self.e, name)
            g = earth(t).observe(planet).astrometric()

            self.eq(ra * tau / 24.0, g.ra, 0.001 * arcsecond)
            self.eq(dec * tau / 360.0, g.dec, 0.001 * arcsecond)
            self.eq(dis, g.distance, 0.1 * meter)

    def test_app_planet(self):

        for t, name in product((T0, TA, TB), planets_to_test):
            obj = c.make_object(0, planet_codes[name], b'planet', None)
            ra, dec, dis = c.app_planet(t, obj)

            earth = self.e.earth
            planet = getattr(self.e, name)
            g = earth(t).observe(planet).apparent()

            self.eq(ra * tau / 24.0, g.ra, 0.001 * arcsecond)
            self.eq(dec * tau / 360.0, g.dec, 0.001 * arcsecond)
            self.eq(dis, g.distance, 0.1 * meter)

    def test_topo_planet(self):
        position = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        ggr = coordinates.Topos('75 W', '45 N', 0.0,
                                temperature=10.0, pressure=1010.0)
        ggr.earth = self.e.earth
        ggr.ephemeris = self.e
        delta_t = 0

        for t, name in product((T0, TA, TB), planets_to_test):
            obj = c.make_object(0, planet_codes[name], b'planet', None)
            ra, dec, dis = c.topo_planet(t, delta_t, obj, position)

            planet = getattr(self.e, name)
            g = ggr(t).observe(planet).apparent()

            self.eq(ra * tau / 24.0, g.ra, 0.001 * arcsecond)
            self.eq(dec * tau / 360.0, g.dec, 0.001 * arcsecond)
            self.eq(dis, g.distance, 0.1 * meter)  # TODO: improve this?

    # Tests of basic functions (in alphabetical order by NOVAS name).

    def test_era(self):
        self.delta = 1e-12
        self.eq(c.era(T0), timescales.earth_rotation_angle(T0))
        self.eq(c.era(TA), timescales.earth_rotation_angle(TA))
        self.eq(c.era(TB), timescales.earth_rotation_angle(TB))

    def test_earth_tilt(self):
        self.delta = 1e-14
        for a, b in zip(c.e_tilt(T0), nutationlib.earth_tilt(A0)):
            self.eq(a, b)
        for a, b in zip(c.e_tilt(TA), nutationlib.earth_tilt(AA)):
            self.eq(a, b)
        for a, b in zip(c.e_tilt(TB), nutationlib.earth_tilt(AB)):
            self.eq(a, b)

    def test_equation_of_the_equinoxes_complimentary_terms(self):
        self.delta = 1e-23

        self.eq(nutationlib.equation_of_the_equinoxes_complimentary_terms(A0),
                c.ee_ct(T0, 0.0, 0))
        self.eq(nutationlib.equation_of_the_equinoxes_complimentary_terms(AA),
                c.ee_ct(TA, 0.0, 0))
        self.eq(nutationlib.equation_of_the_equinoxes_complimentary_terms(AB),
                c.ee_ct(TB, 0.0, 0))

    def test_frame_tie(self):
        self.delta = 1e-15
        v = array((1, 2, 3))

        for a, b in zip(c.frame_tie(v, 0),
                        v.dot(framelib.ICRS_to_J2000)):
            self.eq(a, b)

        for a, b in zip(c.frame_tie(v, -1),
                        v.dot(framelib.J2000_to_ICRS)):
            self.eq(a, b)

    def test_fundamental_arguments(self):
        self.delta = 1e-12

        a = nutationlib.fundamental_arguments(jcentury(T0))
        b = c.fund_args(jcentury(T0))
        for i in range(5):
            self.eq(a[i], b[i])

        a = nutationlib.fundamental_arguments(jcentury(TA))
        b = c.fund_args(jcentury(TA))
        for i in range(5):
            self.eq(a[i], b[i])

        a = nutationlib.fundamental_arguments(jcentury(TB))
        b = c.fund_args(jcentury(TB))
        for i in range(5):
            self.eq(a[i], b[i])

    def test_geo_posvel(self):
        self.delta = 1e-13

        obs1 = c.make_observer_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)
        ggr = coordinates.Topos('75 W', '45 N', 0.0,
                                temperature=10.0, pressure=1010.0)
        delta_t = 0.0

        for v1, v2 in zip(c.geo_posvel(T0, delta_t, obs1),
                          earthlib.geocentric_position_and_velocity(ggr, A0)):
            for a, b in zip(v1, v2):
                self.eq(a, b)

    def test_iau2000a(self):
        self.delta = 1e-19

        self.eq(nutationlib.iau2000a(A0)[0], c.nutation.iau2000a(T0, 0.0)[0])
        self.eq(nutationlib.iau2000a(A0)[1], c.nutation.iau2000a(T0, 0.0)[1])

        self.eq(nutationlib.iau2000a(AA)[0], c.nutation.iau2000a(TA, 0.0)[0])
        self.eq(nutationlib.iau2000a(AA)[1], c.nutation.iau2000a(TA, 0.0)[1])

        self.eq(nutationlib.iau2000a(AB)[0], c.nutation.iau2000a(TB, 0.0)[0])
        self.eq(nutationlib.iau2000a(AB)[1], c.nutation.iau2000a(TB, 0.0)[1])

    def test_mean_obliq(self):
        self.delta = 0

        self.eq(c.mean_obliq(T0), nutationlib.mean_obliquity(T0))
        self.eq(c.mean_obliq(TA), nutationlib.mean_obliquity(TA))
        self.eq(c.mean_obliq(TB), nutationlib.mean_obliquity(TB))

    def test_nutation(self):
        self.delta = 1e-15
        v = array([1, 2, 3])

        for a, b in zip(c_nutation(T0, v, direction=0),
                        v.dot(nutationlib.nutation_matrix(A0)[:,:,0])):
            self.eq(a, b)

        for a, b in zip(c_nutation(TA, v, direction=0),
                        v.dot(nutationlib.nutation_matrix(AA)[:,:,0])):
            self.eq(a, b)

        for a, b in zip(c_nutation(TB, v, direction=1),
                        v.dot(nutationlib.nutation_matrix(AB).T[0])):
            self.eq(a, b)

    def test_precession(self):
        self.delta = 1e-15
        v = array([1, 2, 3])
        c.precession(T0, v, TA)

        for a, b in zip(c.precession(T0, v, TA),
                        v.dot(precessionlib.precession_matrix(TA))):
            self.eq(a, b)

        for a, b in zip(c.precession(TB, v, T0),
                        v.dot(precessionlib.precession_matrix(TB).T)):
            self.eq(a, b)

    def test_sidereal_time(self):
        delta_t = 0.0
        self.delta = 1e-13
        self.eq(c.sidereal_time(T0, 0.0, delta_t, False, True),
                timescales.sidereal_time(T0, delta_t))
        self.eq(c.sidereal_time(TA, 0.0, delta_t, False, True),
                timescales.sidereal_time(TA, delta_t))
        self.eq(c.sidereal_time(TB, 0.0, delta_t, False, True),
                timescales.sidereal_time(TB, delta_t))

    def test_terra(self):
        self.delta = 1e-18

        obs1 = c.make_on_surface(45.0, -75.0, 0.0, 10.0, 1010.0)

        class Obs(object):
            latitude = 45.0 * angles.DEG2RAD
            longitude = -75.0 * angles.DEG2RAD
            elevation = 0.0
        obs2 = Obs()

        for v1, v2 in zip(c.terra(obs1, 11.0), earthlib.terra(obs2, 11.0)):
            for a, b in zip(v1, v2):
                self.eq(a, b)

        for v1, v2 in zip(c.terra(obs1, 23.9), earthlib.terra(obs2, 23.9)):
            for a, b in zip(v1, v2):
                self.eq(a, b)

    def test_tdb2tt(self):
        self.delta = 1e-16
        self.eq(c.tdb2tt(T0)[1], timescales.tdb_minus_tt(T0))
        self.eq(c.tdb2tt(TA)[1], timescales.tdb_minus_tt(TA))
        self.eq(c.tdb2tt(TB)[1], timescales.tdb_minus_tt(TB))

def jcentury(t):
    return (t - T0) / 36525.0
