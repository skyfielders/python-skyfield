"""Rebuild the test data that compares Skyfield to the NOVAS library."""

from __future__ import print_function
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
else:
    from novas.compat import eph_manager
    from novas.constants import T0
    nutation_function = novas.nutation
    # from novas.compat import nutation as nutation_module

def main():
    jd_start, jd_end, number = eph_manager.ephem_open()
    output({}, """\

        import pytest
        from numpy import abs
        from skyfield.api import earth, mars
        from skyfield.jpllib import Ephemeris

        try:
            import de405
            de405 = Ephemeris(de405)
        except ImportError:
            pytestmark = pytest.mark.skipif(True, reason='de405 unavailable')

        arcsecond = 1.0 / 60.0 / 60.0
        ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0

        def compare(value, benchmark_value, tolerance):
            assert abs(value - benchmark_value) < tolerance

        """)
    planet = 'mars'
    code = 4

    moon_landing = novas.julian_date(1969, 7, 20, 20.0 + 18.0/60.0)
    first_hubble_image = novas.julian_date(1990, 5, 20)
    voyager_intersellar = novas.julian_date(2012, 8, 25)

    dates = [moon_landing, first_hubble_image, T0, voyager_intersellar]

    for i, jd in enumerate(dates):
        obj = novas.make_object(0, code, 'planet{}'.format(code), None)
        ra, dec, distance = novas.astro_planet(jd, obj)
        output(locals(), """\

        def test_{planet}{i}():
            a = de405.earth(tt={jd!r}).observe(de405.{planet})
            ra, dec, distance = a.radec()
            compare(ra.hours, {ra!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec!r}, 0.001 * arcsecond)

        """)

def output(dictionary, template):
    print(dedent(template).format(**dictionary).strip('\n'))
    print()

if __name__ == '__main__':
    main()
