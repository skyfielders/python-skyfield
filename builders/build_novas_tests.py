"""Rebuild the test data that compares Skyfield to the NOVAS library."""

from __future__ import print_function
from itertools import product
from sys import exit
from textwrap import dedent

try:
    from novas import compat as novas
    import novas_de405
except ImportError:
    print(dedent("""\
        Error: to rebuild NOVAS test data, you must install both the "novas"
        package and its default ephemeris:

        pip install novas novas_de405

        """))
    exit(2)

from novas.compat import eph_manager
from novas.constants import T0
nutation_function = novas.nutation
# from novas.compat import nutation as nutation_module

planets = [('mercury', 1), ('venus', 2), ('mars', 4), ('jupiter', 5),
           ('saturn', 6), ('uranus', 7), ('neptune', 8), ('pluto', 9),
           ('sun', 10), ('moon', 11)]

def main():
    jd_start, jd_end, number = eph_manager.ephem_open()
    output({}, """\

        import pytest
        from numpy import abs, array, max
        from skyfield import earthlib, nutationlib, timelib
        from skyfield.api import JulianDate, earth, mars
        from skyfield.constants import AU_M
        from skyfield.functions import length_of
        from skyfield.jpllib import Ephemeris


        try:
            import de405
            de405 = Ephemeris(de405)
        except ImportError:
            pytestmark = pytest.mark.skipif(True, reason='de405 unavailable')

        one_second = 1.0 / 24.0 / 60.0 / 60.0
        arcsecond = 1.0 / 60.0 / 60.0
        ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
        meter = 1.0 / AU_M

        def compare(value, expected_value, epsilon):
            if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
                assert max(abs(value - expected_value)) <= epsilon
            else:
                assert abs(value - expected_value) <= epsilon

        """)

    moon_landing = novas.julian_date(1969, 7, 20, 20.0 + 18.0/60.0)
    first_hubble_image = novas.julian_date(1990, 5, 20)
    voyager_intersellar = novas.julian_date(2012, 8, 25)

    date_vector = [moon_landing, first_hubble_image, T0, voyager_intersellar]
    dates = date_vector + [date_vector]

    output_subroutine_tests(dates)
    output_geocentric_tests(dates)
    output_topocentric_tests(dates)


def output_subroutine_tests(dates):
    boring_dates = [d for d in dates if not isinstance(d, list)]

    def shorter_cal_date(jd):
        y, m, d, h = novas.cal_date(jd)
        return y, m, d + h / 24.0 - 0.5

    for i, jd in enumerate(boring_dates):
        cal_date = call(shorter_cal_date, jd)
        output(locals(), """\
            def test_calendar_date_{i}():
                compare(timelib.calendar_date({jd!r}), array({cal_date}), 0.0)
            """)

    for i, jd in enumerate(boring_dates):
        angle = novas.era(jd)
        output(locals(), """\
            def test_earth_rotation_angle_date{i}():
                compare(earthlib.earth_rotation_angle({jd!r}), {angle},
                        0.000001 * arcsecond)
            """)

    for i, jd in enumerate(boring_dates):
        angles = novas.e_tilt(jd)
        output(locals(), """\
            def test_earth_tilt_date{i}():
                compare(nutationlib.earth_tilt(JulianDate(tdb={jd!r})),
                        array({angles}), 0.00001 * arcsecond)
            """)

    for i, jd in enumerate(boring_dates):
        terms = novas.ee_ct(jd, 0.0, 0)
        output(locals(), """\
            def test_equation_of_the_equinoxes_complimentary_terms_date{i}():
                compare(nutationlib.equation_of_the_equinoxes_complimentary_terms({jd!r}),
                        array({terms}), 0.0000000000000001 * arcsecond)
            """)


def output_geocentric_tests(dates):
    for (planet, code), (i, jd) in product(planets, enumerate(dates)):
        obj = novas.make_object(0, code, 'planet{}'.format(code), None)

        ra1, dec1, distance1 = call(novas.astro_planet, jd, obj)
        ra2, dec2, distance2 = call(novas.virtual_planet, jd, obj)
        ra3, dec3, distance3 = call(novas.app_planet, jd, obj)

        assert distance1 == distance2 == distance3

        output(locals(), """\

        def test_{planet}_geocentric_date{i}():
            jd = JulianDate(tt={jd!r})
            e = de405.earth(jd)

            distance = length_of((e - de405.{planet}(jd)).position.AU)
            compare(distance, {distance1!r}, 0.5 * meter)

            astrometric = e.observe(de405.{planet})
            ra, dec, distance = astrometric.radec()
            compare(ra.hours, {ra1!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec1!r}, 0.001 * arcsecond)

            apparent = astrometric.apparent()
            ra, dec, distance = apparent.radec()
            compare(ra.hours, {ra2!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec2!r}, 0.001 * arcsecond)

            ra, dec, distance = apparent.radec(epoch='date')
            compare(ra.hours, {ra3!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec3!r}, 0.001 * arcsecond)

        """)


def output_topocentric_tests(dates):
    usno = novas.make_on_surface(38.9215, -77.0669, 92.0, 10.0, 1010.0)
    for (planet, code), (i, jd) in product(planets, enumerate(dates)):
        obj = novas.make_object(0, code, 'planet{}'.format(code), None)

        ra1, dec1, distance1 = call(novas.local_planet, jd, 0.0, obj, usno)
        ra2, dec2, distance2 = call(novas.topo_planet, jd, 0.0, obj, usno)
        alt, az = call(altaz_maneuver, jd, obj, usno)

        output(locals(), """\

        def test_{planet}_topocentric_date{i}():
            jd = JulianDate(tt={jd!r})
            usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

            apparent = usno(jd).observe(de405.{planet}).apparent()
            ra, dec, distance = apparent.radec()
            compare(ra.hours, {ra1!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec1!r}, 0.001 * arcsecond)

            ra, dec, distance = apparent.radec(epoch='date')
            compare(ra.hours, {ra2!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec2!r}, 0.001 * arcsecond)

            alt, az, distance = apparent.altaz()
            compare(alt.degrees, {alt!r}, 0.001 * arcsecond)
            compare(az.degrees, {az!r}, 0.001 * arcsecond)

        """)


def altaz_maneuver(jd, obj, place):
    """Simplify a pair of complicated USNO calls to a single callable."""
    xp = yp = 0.0
    ra, dec, distance = novas.topo_planet(jd, 0.0, obj, place)
    (zd, az), (ra, dec) = novas.equ2hor(jd, 0.0, xp, yp, place, ra, dec)
    return 90.0 - zd, az


def call(function, jd, *args):
    """Call function either once, or as many times as `jd` dictates."""

    if isinstance(jd, float):
        return function(jd, *args)

    answers = [function(n, *args) for n in jd]
    return zip(*answers)


def output(dictionary, template):
    print(dedent(template).format(**dictionary).strip('\n'))
    print()


if __name__ == '__main__':
    main()
