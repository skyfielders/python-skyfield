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
        from numpy import abs
        from skyfield.api import JulianDate, earth, mars
        from skyfield.constants import AU_M
        from skyfield.functions import length_of
        from skyfield.jpllib import Ephemeris

        try:
            import de405
            de405 = Ephemeris(de405)
        except ImportError:
            pytestmark = pytest.mark.skipif(True, reason='de405 unavailable')

        arcsecond = 1.0 / 60.0 / 60.0
        ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
        meter = 1.0 / AU_M

        def compare(value, benchmark_value, tolerance):
            assert abs(value - benchmark_value) < tolerance

        """)

    moon_landing = novas.julian_date(1969, 7, 20, 20.0 + 18.0/60.0)
    first_hubble_image = novas.julian_date(1990, 5, 20)
    voyager_intersellar = novas.julian_date(2012, 8, 25)

    dates = [moon_landing, first_hubble_image, T0, voyager_intersellar]

    output_geocentric_tests(dates)
    output_topocentric_tests(dates)


def output_geocentric_tests(dates):
    for (planet, code), (i, jd) in product(planets, enumerate(dates)):
        obj = novas.make_object(0, code, 'planet{}'.format(code), None)

        ra1, dec1, distance1 = novas.astro_planet(jd, obj)
        ra2, dec2, distance2 = novas.virtual_planet(jd, obj)
        ra3, dec3, distance3 = novas.app_planet(jd, obj)

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
        xp = yp = 0.0

        ra1, dec1, distance1 = novas.local_planet(jd, 0.0, obj, usno)
        ra2, dec2, distance2 = novas.topo_planet(jd, 0.0, obj, usno)
        (zd, az), (ra, dec) = novas.equ2hor(jd, 0.0, xp, yp, usno, ra2, dec2)
        alt = 90.0 - zd

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


def output(dictionary, template):
    print(dedent(template).format(**dictionary).strip('\n'))
    print()


if __name__ == '__main__':
    main()
