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
from novas.constants import ASEC2RAD, T0
nutation_function = novas.nutation
import novas.compat.nutation as nutation_module

one_second = 1.0 / 24.0 / 60.0 / 60.0
planets = [('mercury', 1), ('venus', 2), ('mars', 4),
           ('jupiter barycenter', 5), ('saturn barycenter', 6),
           ('uranus barycenter', 7), ('neptune barycenter', 8),
           ('pluto barycenter', 9), ('sun', 10), ('moon', 11)]

def main():
    jd_start, jd_end, number = eph_manager.ephem_open()
    output({}, """\
        'Auto-generated accuracy tests vs NOVAS (see build_novas_tests.py).'

        from numpy import abs, array, einsum, max
        from skyfield import (earthlib, framelib, nutationlib, positionlib,
                              precessionlib, starlib, timelib)
        from skyfield.api import JulianDate, Timescale, load
        from skyfield.constants import AU_KM, AU_M
        from skyfield.data import hipparcos
        from skyfield.functions import length_of

        OLD_AU_KM = 149597870.691  # TODO: load from de405
        OLD_AU = AU_KM / OLD_AU_KM

        one_second = 1.0 / 24.0 / 60.0 / 60.0
        arcsecond = 1.0 / 60.0 / 60.0
        ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
        meter = 1.0 / AU_M

        def ts():
            yield Timescale()

        def compare(value, expected_value, epsilon):
            if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
                assert max(abs(value - expected_value)) <= epsilon
            else:
                assert abs(value - expected_value) <= epsilon

        def de405():
            yield load('de405.bsp')

        def earth():
            eph = load('de405.bsp')
            yield eph[399]

        """)

    moon_landing = novas.julian_date(1969, 7, 20, 20.0 + 18.0/60.0)
    first_hubble_image = novas.julian_date(1990, 5, 20)
    voyager_intersellar = novas.julian_date(2012, 8, 25)

    date_vector = [moon_landing, first_hubble_image, T0, voyager_intersellar]
    dates = date_vector + [date_vector]

    output_subroutine_tests(dates)
    output_geocentric_tests(dates)
    output_topocentric_tests(dates)
    output_catalog_tests(dates)


def output_subroutine_tests(dates):
    date_floats = [d for d in dates if not isinstance(d, list)]
    delta_t_floats = [+39.707, +57.1136, +63.8285, +66.7846]

    def shorter_cal_date(jd):
        y, m, d, h = novas.cal_date(jd)
        return y, m, d + h / 24.0 - 0.5

    for i, jd in enumerate(date_floats):
        cal_date = call(shorter_cal_date, jd)
        output(locals(), """\
            def test_calendar_date_{i}():
                compare(timelib.calendar_date({jd!r}), array({cal_date}), 0.0)
            """)

    for i, jd in enumerate(date_floats):
        angle = novas.era(jd)
        output(locals(), """\
            def test_earth_rotation_angle_date{i}():
                compare(earthlib.earth_rotation_angle({jd!r}) * 360.0, {angle!r},
                        0.000001 * arcsecond)
            """)

    for i, jd in enumerate(date_floats):
        angles = novas.e_tilt(jd)
        output(locals(), """\
            def test_earth_tilt_date{i}(ts):
                compare(nutationlib.earth_tilt(ts.tdb({jd!r})),
                        array({angles}), 0.00001 * arcsecond)
            """)

    for i, jd in enumerate(date_floats):
        terms = novas.ee_ct(jd, 0.0, 0)
        output(locals(), """\
            def test_equation_of_the_equinoxes_complimentary_terms_date{i}():
                compare(nutationlib.equation_of_the_equinoxes_complimentary_terms({jd!r}),
                        array({terms!r}), 0.0000000000000001 * arcsecond)
            """)

    vector = (1.1, 1.2, 1.3)
    tie1 = novas.frame_tie(vector, 0)
    tie2 = novas.frame_tie(vector, -1)
    output(locals(), """\
        def test_forward_frame_tie():
            compare(framelib.ICRS_to_J2000.dot({vector}), {tie1}, 1e-15)

        def test_reverse_frame_tie():
            compare(framelib.ICRS_to_J2000.T.dot({vector}), {tie2}, 1e-15)
        """)

    for i, jd in enumerate(date_floats):
        jcentury = (jd - T0) / 36525.0
        arguments = novas.fund_args(jcentury)
        output(locals(), """\
            def test_fundamental_arguments_date{i}():
                compare(nutationlib.fundamental_arguments({jcentury!r}),
                        array({arguments}), 0.000000002 * arcsecond)
            """)

    for i, jd in enumerate(date_floats):
        psi, eps = nutation_module.iau2000a(jd, 0.0)
        psi *= 1e7 / ASEC2RAD
        eps *= 1e7 / ASEC2RAD
        output(locals(), """\
            def test_iau2000a_date{i}():
                compare(nutationlib.iau2000a({jd!r}),
                        array([{psi!r}, {eps!r}]), 0.001)
            """)

    for i, args in enumerate([
          (-4712, 1, 1, 0.0),
          (-4712, 3, 1, 0.0),
          (-4712, 12, 31, 0.5),
          (-241, 3, 25, 19.0),
          (530, 9, 27, 23.5),
          (1976, 3, 7, 12.5),
          (2000, 1, 1, 0.0),
          ]):
        jd = novas.julian_date(*args)
        output(locals(), """\
            def test_julian_date_function_date{i}():
                compare(timelib.julian_date{args}, {jd!r}, 0.0)
            """)

    for i, jd in enumerate(date_floats):
        angle = novas.mean_obliq(jd)
        output(locals(), """\
            def test_mean_obliquity_date{i}():
                compare(nutationlib.mean_obliquity({jd!r}),
                        {angle!r}, 0.0)  # arcseconds
            """)

    for i, jd in enumerate(date_floats):
        vector = [1.1, 1.2, 1.3]
        result = nutation_function(jd, vector)
        output(locals(), """\
            def test_nutation_date{i}(ts):
                matrix = nutationlib.compute_nutation(ts.tdb({jd!r}))
                result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
                compare({result},
                        result, 1e-14)
            """)

    for i, jd in enumerate(date_floats):
        vector = [1.1, 1.2, 1.3]
        result = novas.precession(T0, vector, jd)
        output(locals(), """\
            def test_precession_date{i}():
                matrix = precessionlib.compute_precession({jd!r})
                result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
                compare({result},
                        result, 1e-15)
            """)

    for i, jd in enumerate(date_floats):
        result1 = novas.sidereal_time(jd, 0.0, 0.0, False, True)
        result2 = novas.sidereal_time(jd, 0.0, 99.9, False, True)
        output(locals(), """\
            def test_sidereal_time_on_date{i}():
                jd = Timescale(delta_t=0.0).tt({jd!r})
                compare(earthlib.sidereal_time(jd), {result1!r}, 1e-13)

            def test_sidereal_time_with_nonzero_delta_t_on_date{i}():
                jd = Timescale(delta_t=99.9).tt({jd!r} + 99.9 * one_second)
                compare(earthlib.sidereal_time(jd), {result2!r}, 1e-13)
            """)

    p, v = novas.starvectors(novas.make_cat_entry(
        'POLARIS', 'HIP', 0, 2.530301028, 89.264109444,
        44.22, -11.75, 7.56, -17.4))
    output(locals(), """\
        def test_star_vector():
            star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                                ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                                parallax_mas=7.56, radial_km_per_s=-17.4)
            star.au_km = OLD_AU_KM
            star._compute_vectors()
            compare(star._position_au,
                    {p!r},
                    1e3 * meter)
            compare(star._velocity_au_per_d,
                    {v!r},
                    1e-3 * meter)  # TODO: was 1e-6 before switch to modern au
        """)

    atp = product([-5, -1, 15, 89.95], [10, 25], [1010, 1013.25])

    for i, (angle, temperature, pressure) in enumerate(atp):
        location = novas.make_on_surface(0.0, 0.0, 0, temperature, pressure)
        r = novas.refract(location, 90 - angle, 2)
        output(locals(), """\
            def test_refraction{i}():
                r = earthlib.refraction({angle}, {temperature}, {pressure})
                compare(r, {r!r}, 0.001 * arcsecond)
            """)

    northpole = novas.make_on_surface(90.0, 0.0, 0.0, 10.0, 1010.0)

    for i, angle in enumerate([-90, -2, -1, 0, 1, 3, 9, 90]):
        alt, az = altaz_maneuver(T0, northpole, 0.0, angle, ref=2)
        output(locals(), """\
            def test_refract{i}():
                alt = earthlib.refract({angle!r}, 10.0, 1010.0)
                compare(alt, {alt!r}, 0.000000001 * arcsecond)
            """)

    usno = novas.make_on_surface(38.9215, -77.0669, 92.0, 10.0, 1010.0)

    ra = 12.34
    for i, (tt, dec) in enumerate(product(date_floats, [56.78, -67.89])):
        alt, az = altaz_maneuver(tt, usno, ra, dec, ref=0)
        output(locals(), """\
            def test_from_altaz_{i}(earth):
                jd = Timescale(delta_t=0.0).tt({tt!r})
                usno = earth.topos(
                    '38.9215 N', '77.0669 W', elevation_m=92.0)
                a = usno.at(jd).from_altaz(alt_degrees={alt!r}, az_degrees={az!r})
                ra, dec, distance = a.radec(epoch=jd)
                compare(ra.hours, {ra!r}, 0.000000001 * arcsecond)
                compare(dec.degrees, {dec!r}, 0.000000001 * arcsecond)
            """)

    for i, (tt, delta_t) in enumerate(zip(date_floats, delta_t_floats)):
        jd_low = xp = yp = 0.0
        vector = [1.1, 1.2, 1.3]
        ut1 = tt - delta_t * one_second
        result = novas.ter2cel(ut1, jd_low, delta_t, xp, yp, vector)
        output(locals(), """\
            def test_ITRF_to_GCRS_conversion_on_date{i}():
                jd = Timescale(delta_t={delta_t!r}).tt({tt!r})
                position = positionlib.ITRF_to_GCRS(jd, {vector!r})
                compare(position, {result!r}, 1e-13)
            """)

    for i, jd_tdb in enumerate(date_floats):
        result = novas.tdb2tt(jd_tdb)[1]
        output(locals(), """\
            def test_tdb_minus_tt_on_date{i}():
                result = timelib.tdb_minus_tt({jd_tdb!r})
                compare(result, {result!r}, 1e-16)
            """)


def output_geocentric_tests(dates):
    for (planet, code), (i, jd) in product(planets, enumerate(dates)):
        slug = slugify(planet)
        obj = novas.make_object(0, code, 'planet{}'.format(code), None)

        ra1, dec1, distance1 = call(novas.astro_planet, jd, obj)
        ra2, dec2, distance2 = call(novas.virtual_planet, jd, obj)
        ra3, dec3, distance3 = call(novas.app_planet, jd, obj)

        assert distance1 == distance2 == distance3

        output(locals(), """\

        def test_{slug}_geocentric_date{i}(de405):
            jd = JulianDate(tt={jd!r})
            e = de405['earth'].at(jd)
            p = de405[{planet!r}]

            distance = length_of((e - p.at(jd)).position.au)
            compare(distance * OLD_AU, {distance1!r}, 0.5 * meter)

            astrometric = e.observe(p)
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

    # And, one star.

    polaris = novas.make_cat_entry(
        'POLARIS', 'HIP', 0, 2.530301028, 89.264109444,
        44.22, -11.75, 7.56, -17.4)

    starlist = [('polaris', polaris)]

    for (name, star), (i, jd) in product(starlist, enumerate(dates)):

        ra1, dec1 = call(novas.astro_star, jd, star)
        ra2, dec2 = call(novas.virtual_star, jd, star)
        ra3, dec3 = call(novas.app_star, jd, star)

        output(locals(), """\

        def test_{name}_geocentric_date{i}(earth):
            e = earth.at(Timescale().tt({jd!r}))
            star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                                ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                                parallax_mas=7.56, radial_km_per_s=-17.4)

            astrometric = e.observe(star)
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
        slug = slugify(planet)
        obj = novas.make_object(0, code, 'planet{}'.format(code), None)

        ra1, dec1, distance1 = call(novas.local_planet, jd, 0.0, obj, usno)
        ra2, dec2, distance2 = call(novas.topo_planet, jd, 0.0, obj, usno)
        alt, az = call(altaz_maneuver, jd, usno, ra2, dec2)
        alt2, az2 = call(altaz_maneuver, jd, usno, ra2, dec2, 1)
        alt3, az3 = call(altaz_maneuver, jd, usno, ra2, dec2, 2)

        output(locals(), """\

        def test_{slug}_topocentric_date{i}(de405):
            jd = Timescale(delta_t=0.0).tt({jd!r})
            earth = de405['earth']
            usno = earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

            apparent = usno.at(jd).observe(de405[{planet!r}]).apparent()
            ra, dec, distance = apparent.radec()
            compare(ra.hours, {ra1!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec1!r}, 0.001 * arcsecond)

            ra, dec, distance = apparent.radec(epoch='date')
            compare(ra.hours, {ra2!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec2!r}, 0.001 * arcsecond)

            alt, az, distance = apparent.altaz()
            compare(alt.degrees, {alt!r}, 0.001 * arcsecond)
            compare(az.degrees, {az!r}, 0.001 * arcsecond)

            alt, az, distance = apparent.altaz('standard')
            compare(alt.degrees, {alt2!r}, 0.001 * arcsecond)
            compare(az.degrees, {az2!r}, 0.001 * arcsecond)

            alt, az, distance = apparent.altaz(10.0, 1010.0)
            compare(alt.degrees, {alt3!r}, 0.001 * arcsecond)
            compare(az.degrees, {az3!r}, 0.001 * arcsecond)

        """)


def output_catalog_tests(dates):
    polaris1991 = novas.make_cat_entry(
        '11767', 'HIP', 0, 0.0, 89.26413805,
        44.22, -11.74, 7.56, 0.0)
    polaris1991.ra = 37.94614689  # HIP uses degrees, not hours
    polaris = novas.transform_hip(polaris1991)
    for i, jd in enumerate(dates):
        ra, dec = call(novas.astro_star, jd, polaris)
        output(locals(), r"""

        def test_hipparcos_conversion{i}(earth):
            line = 'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
            star = hipparcos.parse(line)
            compare(star.ra.hours, {polaris.ra!r}, 0.001 * ra_arcsecond)
            compare(star.dec.degrees, {polaris.dec!r}, 0.001 * arcsecond)
            ra, dec, distance = earth.at(Timescale().tt({jd})).observe(star).radec()
            compare(ra.hours, {ra!r}, 0.001 * ra_arcsecond)
            compare(dec.degrees, {dec!r}, 0.001 * arcsecond)

        """)


def altaz_maneuver(jd, place, ra, dec, ref=0):
    """Wrapper that simplifies a complicated USNO call."""
    xp = yp = 0.0
    (zd, az), (ra, dec) = novas.equ2hor(jd, 0.0, xp, yp, place, ra, dec, ref)
    return 90.0 - zd, az


def call(function, *args):
    """Call function as many times as any array arguments dictate."""

    length = max(len(arg) if hasattr(arg, '__len__') else 0 for arg in args)
    if not length:
        return function(*args)
    argstacks = [(arg if hasattr(arg, '__len__') else [arg] * length)
                 for arg in args]
    answers = [function(*arglist) for arglist in zip(*argstacks)]
    return list(zip(*answers))


def slugify(name):
    """Turn 'jupiter_barycenter' into 'jupiter barycenter'."""
    return name.replace(' ', '_')


def output(dictionary, template):
    print(dedent(template).format(**dictionary).strip('\n'))
    print()


if __name__ == '__main__':
    main()
