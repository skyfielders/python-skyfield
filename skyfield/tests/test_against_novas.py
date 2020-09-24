'Auto-generated accuracy tests vs NOVAS (see build_novas_tests.py).'

from numpy import abs, array, einsum, max
from skyfield import (earthlib, framelib, nutationlib, positionlib,
                      precessionlib, starlib, timelib)
from skyfield.api import Topos, load
from skyfield.constants import AU_KM, AU_M
from skyfield.data import hipparcos
from skyfield.functions import BytesIO, length_of
from .fixes import low_precision_ERA

OLD_AU_KM = 149597870.691  # TODO: load from de405
OLD_AU = AU_KM / OLD_AU_KM

one_second = 1.0 / 24.0 / 60.0 / 60.0
arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
meter = 1.0 / AU_M

def ts():
    yield load.timescale()

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

def reduce_precision(t):
    # The NOVAS library uses only 64-bit precision for Julian dates.
    # Some of these tests are so sensitive they can see the difference!
    # So we need to collapse the Julian dates back into single floats.
    delta = t.tdb - t.tt
    t.whole = t.tdb
    t.tt_fraction = delta
    t.tdb_fraction = 0.0

def test_calendar_date_0():
    compare(timelib.compute_calendar_date(2440423), array((1969, 7, 20.0)), 0.0)

def test_calendar_date_1():
    compare(timelib.compute_calendar_date(2448031), array((1990, 5, 19.0)), 0.0)

def test_calendar_date_2():
    compare(timelib.compute_calendar_date(2451545), array((2000, 1, 1.0)), 0.0)

def test_calendar_date_3():
    compare(timelib.compute_calendar_date(2456164), array((2012, 8, 24.0)), 0.0)

def test_earth_rotation_angle_date0():
    compare(earthlib.earth_rotation_angle(2440423.345833333) * 360.0, 243.3216078027496,
            0.000001 * arcsecond)

def test_earth_rotation_angle_date1():
    compare(earthlib.earth_rotation_angle(2448031.5) * 360.0, 237.5118441792128,
            0.000001 * arcsecond)

def test_earth_rotation_angle_date2():
    compare(earthlib.earth_rotation_angle(2451545.0) * 360.0, 280.46061837504,
            0.000001 * arcsecond)

def test_earth_rotation_angle_date3():
    compare(earthlib.earth_rotation_angle(2456164.5) * 360.0, 333.4965831957672,
            0.000001 * arcsecond)

def test_earth_tilt_date0(ts):
    compare(nutationlib.earth_tilt(ts.tdb_jd(2440423.345833333)),
            array((23.443240959852666, 23.445702723464045, 0.15929455696954214, 2.604727521416375, 8.862349000962691)), 0.00001 * arcsecond)

def test_earth_tilt_date1(ts):
    compare(nutationlib.earth_tilt(ts.tdb_jd(2448031.5)),
            array((23.440530953006782, 23.442178709915066, 0.7110982205507752, 11.628148141964171, 5.931924869819427)), 0.00001 * arcsecond)

def test_earth_tilt_date2(ts):
    compare(nutationlib.earth_tilt(ts.tdb_jd(2451545.0)),
            array((23.439279444444445, 23.437676833867652, -0.852016747090803, -13.931996330960066, -5.769398076465291)), 0.00001 * arcsecond)

def test_earth_tilt_date3(ts):
    compare(nutationlib.earth_tilt(ts.tdb_jd(2456164.5)),
            array((23.43763397776759, 23.43645066577372, 0.977087608170215, 15.976729533480038, -4.259923177932873)), 0.00001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date0():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2440423.345833333),
            array(-1.4592438843164885e-09), 0.0000000000000001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date1():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2448031.5),
            array(-9.909270679336256e-09), 0.0000000000000001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date2():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2451545.0),
            array(1.021330096302465e-08), 0.0000000000000001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date3():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2456164.5),
            array(-1.082315527387237e-08), 0.0000000000000001 * arcsecond)

def test_forward_frame_tie():
    compare(framelib.ICRS_to_J2000.dot((1.1, 1.2, 1.3)), (1.100000019790573, 1.2000001208396125, 1.2999998717098593), 1e-15)

def test_reverse_frame_tie():
    compare(framelib.ICRS_to_J2000.T.dot((1.1, 1.2, 1.3)), (1.0999999802094143, 1.1999998791603803, 1.300000128290131), 1e-15)

def test_fundamental_arguments_date0():
    compare(nutationlib.fundamental_arguments(-0.3044942961441969),
            array((-1.559784616935014, -2.8619278194907483, -2.7748368269156427, -4.947060102171707, 6.178085194718492)), 0.000000002 * arcsecond)

def test_fundamental_arguments_date1():
    compare(nutationlib.fundamental_arguments(-0.09619438740588637),
            array((-0.8532784044768771, -3.933579124091533, -5.376486844354831, -0.9485312704748627, 5.429677887938805)), 0.000000002 * arcsecond)

def test_fundamental_arguments_date2():
    compare(nutationlib.fundamental_arguments(0.0),
            array((2.355555743493879, 6.24006012692298, 1.6279050815375191, 5.198466588650503, 2.182439196615671)), 0.000000002 * arcsecond)

def test_fundamental_arguments_date3():
    compare(nutationlib.fundamental_arguments(0.12647501711156742),
            array((0.15181719486225662, 4.023151622222436, 0.10917837795937814, 1.6234303368860354, -2.086983188457769)), 0.000000002 * arcsecond)

def test_iau2000a_date0():
    compare(nutationlib.iau2000a(2440423.345833333),
            array([26047275.214163747, 88623490.00962691]), 0.001)

def test_iau2000a_date1():
    compare(nutationlib.iau2000a(2448031.5),
            array([116281481.4196417, 59319248.69819427]), 0.001)

def test_iau2000a_date2():
    compare(nutationlib.iau2000a(2451545.0),
            array([-139319963.30960065, -57693980.764652915]), 0.001)

def test_iau2000a_date3():
    compare(nutationlib.iau2000a(2456164.5),
            array([159767295.3348004, -42599231.779328726]), 0.001)

def test_iau2000b_date0():
    compare(nutationlib.iau2000b(2440423.345833333),
            array([26048264.528388523, 88619675.68529966]), 0.001)

def test_iau2000b_date1():
    compare(nutationlib.iau2000b(2448031.5),
            array([116274598.48837188, 59322174.624764845]), 0.001)

def test_iau2000b_date2():
    compare(nutationlib.iau2000b(2451545.0),
            array([-139316638.88969785, -57694170.77292847]), 0.001)

def test_iau2000b_date3():
    compare(nutationlib.iau2000b(2456164.5),
            array([159765584.29895684, -42598702.03944705]), 0.001)

def test_julian_date_function_date0():
    compare(timelib.julian_date(-4712, 1, 1, 0.0), 37.5, 0.0)

def test_julian_date_function_date1():
    compare(timelib.julian_date(-4712, 3, 1, 0.0), 97.5, 0.0)

def test_julian_date_function_date2():
    compare(timelib.julian_date(-4712, 12, 31, 0.5), 402.5208333333333, 0.0)

def test_julian_date_function_date3():
    compare(timelib.julian_date(-241, 3, 25, 19.0), 1633120.2916666667, 0.0)

def test_julian_date_function_date4():
    compare(timelib.julian_date(530, 9, 27, 23.5), 1914908.4791666667, 0.0)

def test_julian_date_function_date5():
    compare(timelib.julian_date(1976, 3, 7, 12.5), 2442845.0208333335, 0.0)

def test_julian_date_function_date6():
    compare(timelib.julian_date(2000, 1, 1, 0.0), 2451544.5, 0.0)

def test_mean_obliquity_date0():
    compare(nutationlib.mean_obliquity(2440423.345833333),
            84395.6674554696, 0.0)  # arcseconds

def test_mean_obliquity_date1():
    compare(nutationlib.mean_obliquity(2448031.5),
            84385.91143082442, 0.0)  # arcseconds

def test_mean_obliquity_date2():
    compare(nutationlib.mean_obliquity(2451545.0),
            84381.406, 0.0)  # arcseconds

def test_mean_obliquity_date3():
    compare(nutationlib.mean_obliquity(2456164.5),
            84375.48231996332, 0.0)  # arcseconds

def test_nutation_date0(ts):
    matrix = nutationlib.compute_nutation(ts.tdb_jd(2440423.345833333))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0999795659425045, 1.1999568871469584, 1.3000570847072532),
            result, 1e-14)

def test_nutation_date1(ts):
    matrix = nutationlib.compute_nutation(ts.tdb_jd(2448031.5))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0999087778623433, 1.2000195046911977, 1.300059178938428),
            result, 1e-14)

def test_nutation_date2(ts):
    matrix = nutationlib.compute_nutation(ts.tdb_jd(2451545.0))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.1001092900321017, 1.1999681897164485, 1.2999368806421698),
            result, 1e-14)

def test_nutation_date3(ts):
    matrix = nutationlib.compute_nutation(ts.tdb_jd(2456164.5))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0998746654010052, 1.2001050177909849, 1.3000091025381042),
            result, 1e-14)

def test_precession_date0():
    matrix = precessionlib.compute_precession(2440423.345833333)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.1119856573552391, 1.1924703352076302, 1.296727572578774),
            result, 1e-15)

def test_precession_date1():
    matrix = precessionlib.compute_precession(2448031.5)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.1037931405410017, 1.1976299348492718, 1.2989700697273823),
            result, 1e-15)

def test_precession_date2():
    matrix = precessionlib.compute_precession(2451545.0)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.1, 1.1999999999999997, 1.2999999999999998),
            result, 1e-15)

def test_precession_date3():
    matrix = precessionlib.compute_precession(2456164.5)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0950034772583117, 1.203103909268923, 1.3013486728367767),
            result, 1e-15)

def test_sidereal_time_on_date0():
    jd = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    compare(earthlib.sidereal_time(jd), 16.195436227057314, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date0():
    jd = load.timescale(delta_t=99.9).tt_jd(2440423.345833333 + 99.9 * one_second)
    with low_precision_ERA():
        compare(earthlib.sidereal_time(jd), 16.195436229760602, 1e-13)

def test_sidereal_time_on_date1():
    jd = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    compare(earthlib.sidereal_time(jd), 15.825907460288224, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date1():
    jd = load.timescale(delta_t=99.9).tt_jd(2448031.5 + 99.9 * one_second)
    with low_precision_ERA():
        compare(earthlib.sidereal_time(jd), 15.825907462991848, 1e-13)

def test_sidereal_time_on_date2():
    jd = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    compare(earthlib.sidereal_time(jd), 18.69737482696563, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date2():
    jd = load.timescale(delta_t=99.9).tt_jd(2451545.0 + 99.9 * one_second)
    with low_precision_ERA():
        compare(earthlib.sidereal_time(jd), 18.69737482966941, 1e-13)

def test_sidereal_time_on_date3():
    jd = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    compare(earthlib.sidereal_time(jd), 22.243908497165812, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date3():
    jd = load.timescale(delta_t=99.9).tt_jd(2456164.5 + 99.9 * one_second)
    with low_precision_ERA():
        compare(earthlib.sidereal_time(jd), 22.2439084998698, 1e-13)

def test_star_vector():
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)
    star.au_km = OLD_AU_KM
    star._compute_vectors()
    compare(star._position_au,
            (276301.52367964364, 215517.39549460335, 27281454.18783122),
            1e3 * meter)
    compare(star._velocity_au_per_d,
            (-0.006595734315371155, 0.015163885823867606, -0.010102577482634966),
            1e-3 * meter)  # TODO: was 1e-6 before switch to modern au

def test_refraction0():
    r = earthlib.refraction(-5, 10, 1010)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction1():
    r = earthlib.refraction(-5, 10, 1013.25)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction2():
    r = earthlib.refraction(-5, 25, 1010)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction3():
    r = earthlib.refraction(-5, 25, 1013.25)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction4():
    r = earthlib.refraction(-1, 10, 1010)
    compare(r, 0.8296919418249878, 1e-9 * arcsecond)

def test_refraction5():
    r = earthlib.refraction(-1, 10, 1013.25)
    compare(r, 0.8323617426278902, 1e-9 * arcsecond)

def test_refraction6():
    r = earthlib.refraction(-1, 25, 1010)
    compare(r, 0.7879289246190321, 1e-9 * arcsecond)

def test_refraction7():
    r = earthlib.refraction(-1, 25, 1013.25)
    compare(r, 0.7904643394754796, 1e-9 * arcsecond)

def test_refraction8():
    r = earthlib.refraction(15, 10, 1010)
    compare(r, 0.06056215494995108, 1e-9 * arcsecond)

def test_refraction9():
    r = earthlib.refraction(15, 10, 1013.25)
    compare(r, 0.06075703317132469, 1e-9 * arcsecond)

def test_refraction10():
    r = earthlib.refraction(15, 25, 1010)
    compare(r, 0.057513724331664955, 1e-9 * arcsecond)

def test_refraction11():
    r = earthlib.refraction(15, 25, 1013.25)
    compare(r, 0.057698793246593584, 1e-9 * arcsecond)

def test_refraction12():
    r = earthlib.refraction(89.95, 10, 1010)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction13():
    r = earthlib.refraction(89.95, 10, 1013.25)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction14():
    r = earthlib.refraction(89.95, 25, 1010)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refraction15():
    r = earthlib.refraction(89.95, 25, 1013.25)
    compare(r, 0.0, 1e-9 * arcsecond)

def test_refract0():
    alt = earthlib.refract(-90, 10.0, 1010.0)
    compare(alt, -90.0, 1e-9 * arcsecond)

def test_refract1():
    alt = earthlib.refract(-2, 10.0, 1010.0)
    compare(alt, -2.0, 1e-9 * arcsecond)

def test_refract2():
    alt = earthlib.refract(-1, 10.0, 1010.0)
    compare(alt, -0.34540033564054795, 1e-9 * arcsecond)

def test_refract3():
    alt = earthlib.refract(0, 10.0, 1010.0)
    compare(alt, 0.4819388815393779, 1e-9 * arcsecond)

def test_refract4():
    alt = earthlib.refract(1, 10.0, 1010.0)
    compare(alt, 1.362447444478633, 1e-9 * arcsecond)

def test_refract5():
    alt = earthlib.refract(3, 10.0, 1010.0)
    compare(alt, 3.227564692764261, 1e-9 * arcsecond)

def test_refract6():
    alt = earthlib.refract(9, 10.0, 1010.0)
    compare(alt, 9.098059272393698, 1e-9 * arcsecond)

def test_refract7():
    alt = earthlib.refract(90, 10.0, 1010.0)
    compare(alt, 90.0, 1e-9 * arcsecond)

def test_from_altaz_0(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=68.12871390985244, az_degrees=28.979244220884173)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, 56.78, 1e-9 * arcsecond)

def test_from_altaz_1(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=-17.792497521318964, az_degrees=172.51742180816711)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, -67.89, 1e-9 * arcsecond)

def test_from_altaz_2(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=65.8650913573598, az_degrees=34.158756360615946)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, 56.78, 1e-9 * arcsecond)

def test_from_altaz_3(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=-18.43186389552551, az_degrees=170.42969631720953)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, -67.89, 1e-9 * arcsecond)

def test_from_altaz_4(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=68.47898348962792, az_degrees=332.05109419434154)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, 56.78, 1e-9 * arcsecond)

def test_from_altaz_5(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=-17.699091955922242, az_degrees=187.12243108963492)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, -67.89, 1e-9 * arcsecond)

def test_from_altaz_6(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=41.36529829114181, az_degrees=316.19259712235026)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, 56.78, 1e-9 * arcsecond)

def test_from_altaz_7(earth):
    jd = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    usno = earth + Topos(
        '38.9215 N', '77.0669 W', elevation_m=92.0)
    a = usno.at(jd).from_altaz(alt_degrees=-29.282626410822033, az_degrees=204.1557062303077)
    ra, dec, distance = a.radec(epoch=jd)
    compare(ra.hours, 12.34, 1e-9 * arcsecond)
    compare(dec.degrees, -67.89, 1e-9 * arcsecond)

def test_ITRF_to_GCRS_conversion_on_date0():
    jd = load.timescale(delta_t=39.707).tt_jd(2440423.345833333)
    with low_precision_ERA():
        position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (0.5701172053658128, -1.5232987806096392, 1.3017400651201707), 1e-13)

def test_ITRF_to_GCRS_conversion_on_date1():
    jd = load.timescale(delta_t=57.1136).tt_jd(2448031.5)
    with low_precision_ERA():
        position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (0.41362649279562963, -1.5741081933652488, 1.3004216700893525), 1e-13)

def test_ITRF_to_GCRS_conversion_on_date2():
    jd = load.timescale(delta_t=63.8285).tt_jd(2451545.0)
    with low_precision_ERA():
        position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (1.3757008573963405, -0.8702954291925735, 1.3000126987400913), 1e-13)

def test_ITRF_to_GCRS_conversion_on_date3():
    jd = load.timescale(delta_t=66.7846).tt_jd(2456164.5)
    with low_precision_ERA():
        position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (1.5243574049688486, 0.5755748855663746, 1.2980940077752074), 1e-13)

def test_tdb_minus_tt_on_date0():
    result = timelib.tdb_minus_tt(2440423.345833333)
    compare(result, -0.00046798717637519603, 1e-16)

def test_tdb_minus_tt_on_date1():
    result = timelib.tdb_minus_tt(2448031.5)
    compare(result, 0.0011585185926349208, 1e-16)

def test_tdb_minus_tt_on_date2():
    result = timelib.tdb_minus_tt(2451545.0)
    compare(result, -9.575743486095212e-05, 1e-16)

def test_tdb_minus_tt_on_date3():
    result = timelib.tdb_minus_tt(2456164.5)
    compare(result, -0.001241030165936046, 1e-16)

def test_position_and_velocity(de405, ts):
    t = ts.tt_jd(2451545.0)
    e = de405['earth'].at(t)
    compare(e.position.au, (-0.18427156190936703, 0.8847815016994275, 0.3838199425739002), 10 * meter)
    compare(e.velocity.au_per_d, (-0.017202246596776286, -0.0029049259923337044, -0.0012594278596487706), 1e-5 * meter / one_second)

def test_mercury_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mercury']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.3278115470600746, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 7.905384000977572, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 22.332364359841474, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.904987228126012, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 22.333433087908823, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.874971625095716, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 22.415970392044656, 0.0001 * arcsecond)

def test_mercury_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mercury']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.6507044512046538, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.4704717994133576, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 11.2501328449305, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.4701282535729665, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 11.248550502940756, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.4616767226464757, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 11.207785493244957, 0.0001 * arcsecond)

def test_mercury_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mercury']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.4155249674526948, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.13892977357885, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -24.42032494108073, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.13851035907211, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -24.420393338459686, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.138225455402914, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -24.418845803732086, 0.0001 * arcsecond)

def test_mercury_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mercury']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.1264323486728112, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 9.295934662566733, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 16.68579742896488, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 9.295575039086721, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 16.687409731964937, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 9.307566088097714, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 16.631743449679668, 0.0001 * arcsecond)

def test_mercury_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mercury']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (1.3278115470600746, 0.6507044512046538, 1.4155249674526948, 1.1264323486728112), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (7.905384000977572, 2.4704717994133576, 18.13892977357885, 9.295934662566733), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (22.332364359841474, 11.2501328449305, -24.42032494108073, 16.68579742896488), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (7.904987228126012, 2.4701282535729665, 18.13851035907211, 9.295575039086721), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (22.333433087908823, 11.248550502940756, -24.420393338459686, 16.687409731964937), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (7.874971625095716, 2.4616767226464757, 18.138225455402914, 9.307566088097714), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (22.415970392044656, 11.207785493244957, -24.418845803732086, 16.631743449679668), 0.0001 * arcsecond)

def test_venus_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['venus']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.9646045654448725, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 4.966946050917652, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 20.210417323471006, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.966656139420439, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 20.210145917097474, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.93668626355443, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 20.166644671858105, 0.0001 * arcsecond)

def test_venus_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['venus']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.0711674186789975, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 1.161811406279447, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 5.32829157368082, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.1615415906820667, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 5.326768071513868, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.1534174784892788, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 5.277365365528824, 0.0001 * arcsecond)

def test_venus_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['venus']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.1376890757925104, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 15.993350650200568, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -18.451653207795236, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.993038357924485, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -18.450881488018126, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.992790109710333, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -18.44871897642583, 0.0001 * arcsecond)

def test_venus_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['venus']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.7824924286112764, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 7.175585125577371, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 19.874130272238094, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.175312328808404, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 19.87477997549141, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.188033727750362, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 19.85167856390226, 0.0001 * arcsecond)

def test_venus_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['venus']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (0.9646045654448725, 1.0711674186789975, 1.1376890757925104, 0.7824924286112764), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (4.966946050917652, 1.161811406279447, 15.993350650200568, 7.175585125577371), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (20.210417323471006, 5.32829157368082, -18.451653207795236, 19.874130272238094), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (4.966656139420439, 1.1615415906820667, 15.993038357924485, 7.175312328808404), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (20.210145917097474, 5.326768071513868, -18.450881488018126, 19.87477997549141), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (4.93668626355443, 1.1534174784892788, 15.992790109710333, 7.188033727750362), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (20.166644671858105, 5.277365365528824, -18.44871897642583, 19.85167856390226), 0.0001 * arcsecond)

def test_mars_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mars']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.5912188976380217, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 16.0296606272219, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -24.127310308581468, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.02988433983068, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -24.128202621801755, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.99950982315885, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -24.046277103674843, 0.0001 * arcsecond)

def test_mars_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mars']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.430250679602913, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 23.545034875459514, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -4.8822490432210355, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.544892038854186, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -4.88299363089811, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.536847630733252, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -4.935089760397492, 0.0001 * arcsecond)

def test_mars_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mars']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.8496039270835372, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 22.034936616343344, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -13.18070741103498, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.03468760932384, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -13.182134899635477, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.034417492807563, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -13.182689288940116, 0.0001 * arcsecond)

def test_mars_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mars']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.7665523168668773, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 13.894324196598355, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -12.122808318928707, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.894132382683363, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -12.121796956140246, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.9057161859901, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -12.184654273116957, 0.0001 * arcsecond)

def test_mars_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['mars']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (0.5912188976380217, 1.430250679602913, 1.8496039270835372, 1.7665523168668773), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (16.0296606272219, 23.545034875459514, 22.034936616343344, 13.894324196598355), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-24.127310308581468, -4.8822490432210355, -13.18070741103498, -12.122808318928707), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (16.02988433983068, 23.544892038854186, 22.03468760932384, 13.894132382683363), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-24.128202621801755, -4.88299363089811, -13.182134899635477, -12.121796956140246), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.99950982315885, 23.536847630733252, 22.034417492807563, 13.9057161859901), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-24.046277103674843, -4.935089760397492, -13.182689288940116, -12.184654273116957), 0.0001 * arcsecond)

def test_jupiter_barycenter_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['jupiter barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 5.8416003192317465, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.104091505864654, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 0.6513409058207986, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.103936313614676, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 0.6524656208782568, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.07798204538282, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 0.8216129394812305, 0.0001 * arcsecond)

def test_jupiter_barycenter_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['jupiter barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 5.913287883102948, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 6.765154678701348, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 23.170397700122013, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 6.764854244708427, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 23.170736332068763, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 6.755383083025232, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 23.182684693676578, 0.0001 * arcsecond)

def test_jupiter_barycenter_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['jupiter barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 4.621126565890217, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 1.5913207023268698, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 8.595887646396902, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.5914167941833441, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 8.59631203599914, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.5911888424331277, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 8.594250857972387, 0.0001 * arcsecond)

def test_jupiter_barycenter_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['jupiter barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 5.129958529243068, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 4.822841055032964, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 21.649994488649476, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.822764769132736, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 21.64994169521302, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.835670404865468, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 21.67058638943795, 0.0001 * arcsecond)

def test_jupiter_barycenter_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['jupiter barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (5.8416003192317465, 5.913287883102948, 4.621126565890217, 5.129958529243068), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.104091505864654, 6.765154678701348, 1.5913207023268698, 4.822841055032964), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (0.6513409058207986, 23.170397700122013, 8.595887646396902, 21.649994488649476), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.103936313614676, 6.764854244708427, 1.5914167941833441, 4.822764769132736), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (0.6524656208782568, 23.170736332068763, 8.59631203599914, 21.64994169521302), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.07798204538282, 6.755383083025232, 1.5911888424331277, 4.835670404865468), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (0.8216129394812305, 23.182684693676578, 8.594250857972387, 21.67058638943795), 0.0001 * arcsecond)

def test_saturn_barycenter_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['saturn barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 9.382032444401025, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.4627748852420206, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 12.045819985925936, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.462707593703528, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 12.045735497802628, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.4352879582290177, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 11.9115661075769, 0.0001 * arcsecond)

def test_saturn_barycenter_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['saturn barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 9.420484451056101, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 19.814248756112033, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -20.933390198050763, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.81446344451556, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -20.932846451357463, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.805277718955743, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -20.958164640919687, 0.0001 * arcsecond)

def test_saturn_barycenter_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['saturn barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 8.652750126001484, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.584400980536592, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 12.616288735770384, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.584593321351076, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 12.616983167644802, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.584361121508456, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 12.614774672730574, 0.0001 * arcsecond)

def test_saturn_barycenter_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['saturn barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 10.326368974662916, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 13.628484577191722, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -7.659435207931653, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.62827504244793, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -7.658028344724226, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.639628746850631, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -7.723201642102626, 0.0001 * arcsecond)

def test_saturn_barycenter_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['saturn barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (9.382032444401025, 9.420484451056101, 8.652750126001484, 10.326368974662916), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (2.4627748852420206, 19.814248756112033, 2.584400980536592, 13.628484577191722), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (12.045819985925936, -20.933390198050763, 12.616288735770384, -7.659435207931653), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (2.462707593703528, 19.81446344451556, 2.584593321351076, 13.62827504244793), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (12.045735497802628, -20.932846451357463, 12.616983167644802, -7.658028344724226), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (2.4352879582290177, 19.805277718955743, 2.584361121508456, 13.639628746850631), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (11.9115661075769, -20.958164640919687, 12.614774672730574, -7.723201642102626), 0.0001 * arcsecond)

def test_uranus_barycenter_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['uranus barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 18.75197906203834, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.087167068351334, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 0.20723926118363256, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.087010426255667, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 0.20832526777272883, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.061052547705433, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 0.37749969290358576, 0.0001 * arcsecond)

def test_uranus_barycenter_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['uranus barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 18.622417009295177, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.668551452013403, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -23.437331340689163, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.668859170516964, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -23.437016930580615, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.65936113308538, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -23.447681812488984, 0.0001 * arcsecond)

def test_uranus_barycenter_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['uranus barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 20.727159134679393, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 21.165586867541418, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -17.018831731314233, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 21.165269485049027, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -17.020267168405784, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 21.164987614252272, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -17.020320613172004, 0.0001 * arcsecond)

def test_uranus_barycenter_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['uranus barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 19.234768680195387, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 0.4891643148564316, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 2.3565095329111823, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 0.4894463256538988, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 2.358369638516312, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 0.5005500654503398, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 2.429779341040803, 0.0001 * arcsecond)

def test_uranus_barycenter_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['uranus barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (18.75197906203834, 18.622417009295177, 20.727159134679393, 19.234768680195387), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.087167068351334, 18.668551452013403, 21.165586867541418, 0.4891643148564316), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (0.20723926118363256, -23.437331340689163, -17.018831731314233, 2.3565095329111823), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.087010426255667, 18.668859170516964, 21.165269485049027, 0.4894463256538988), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (0.20832526777272883, -23.437016930580615, -17.020267168405784, 2.358369638516312), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.061052547705433, 18.65936113308538, 21.164987614252272, 0.5005500654503398), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (0.37749969290358576, -23.447681812488984, -17.020320613172004, 2.429779341040803), 0.0001 * arcsecond)

def test_neptune_barycenter_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['neptune barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 29.83221264621946, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 15.637210587139663, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -17.67999613660563, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.63739098768298, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -17.68045373026462, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.608486730597075, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -17.583793285519313, 0.0001 * arcsecond)

def test_neptune_barycenter_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['neptune barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 29.490001740438892, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 19.03623522579387, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -21.792864018500975, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.036513633320563, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -21.79251066237039, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.02716408230529, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -21.808047913986808, 0.0001 * arcsecond)

def test_neptune_barycenter_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['neptune barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 31.024491920354496, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 20.362841834121518, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -19.21242523937633, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 20.362475439010588, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -19.213645950878377, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 20.36218815756048, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -19.21323379889766, 0.0001 * arcsecond)

def test_neptune_barycenter_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['neptune barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 28.984118029716345, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 22.252468120719442, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -11.504657215501584, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.252825961036415, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -11.50264948264589, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.2643158309744, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -11.437330191299896, 0.0001 * arcsecond)

def test_neptune_barycenter_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['neptune barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (29.83221264621946, 29.490001740438892, 31.024491920354496, 28.984118029716345), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (15.637210587139663, 19.03623522579387, 20.362841834121518, 22.252468120719442), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-17.67999613660563, -21.792864018500975, -19.21242523937633, -11.504657215501584), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (15.63739098768298, 19.036513633320563, 20.362475439010588, 22.252825961036415), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-17.68045373026462, -21.79251066237039, -19.213645950878377, -11.50264948264589), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.608486730597075, 19.02716408230529, 20.36218815756048, 22.2643158309744), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-17.583793285519313, -21.808047913986808, -19.21323379889766, -11.437330191299896), 0.0001 * arcsecond)

def test_pluto_barycenter_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['pluto barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 32.312971776632494, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.015311208821212, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 16.620557180992588, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.01514128380381, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 16.622990160668607, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 11.989232654068259, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 16.792242650891875, 0.0001 * arcsecond)

def test_pluto_barycenter_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['pluto barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 28.707485955458118, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 15.216302246424346, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -1.3346560528819575, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.216661036271791, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -1.3358630622052712, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.208581663980876, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -1.3022394883151638, 0.0001 * arcsecond)

def test_pluto_barycenter_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['pluto barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 31.064412196006614, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 16.761873062250743, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -11.39643313463007, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.761526675406767, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -11.396301545071504, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.761277438459963, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -11.39428873441123, 0.0001 * arcsecond)

def test_pluto_barycenter_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['pluto barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 31.69909782133193, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.488351288595236, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -19.55219099488885, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.488573622605898, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -19.551729414764313, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.501338273669152, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -19.541227909743732, 0.0001 * arcsecond)

def test_pluto_barycenter_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['pluto barycenter']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (32.312971776632494, 28.707485955458118, 31.064412196006614, 31.69909782133193), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.015311208821212, 15.216302246424346, 16.761873062250743, 18.488351288595236), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (16.620557180992588, -1.3346560528819575, -11.39643313463007, -19.55219099488885), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.01514128380381, 15.216661036271791, 16.761526675406767, 18.488573622605898), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (16.622990160668607, -1.3358630622052712, -11.396301545071504, -19.551729414764313), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (11.989232654068259, 15.208581663980876, 16.761277438459963, 18.501338273669152), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (16.792242650891875, -1.3022394883151638, -11.39428873441123, -19.541227909743732), 0.0001 * arcsecond)

def test_sun_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['sun']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.0160878650466754, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 8.03008088792976, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 20.496475643233936, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 8.02969030304998, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 20.497605463260726, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 8.000108116572395, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 20.58493093599605, 0.0001 * arcsecond)

def test_sun_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['sun']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.0118605934887042, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 3.776110727862678, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 19.907832379364574, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 3.775721385487214, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 19.906601181542, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 3.7666292045824337, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 19.879173772309745, 0.0001 * arcsecond)

def test_sun_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['sun']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.9833276788862821, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.752544254682526, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -23.033309607967187, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.752126228091367, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -23.03376015263556, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.75183797477899, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -23.032488638722818, 0.0001 * arcsecond)

def test_sun_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['sun']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 1.0107820040799866, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 10.268162490439073, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 10.751933902906119, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 10.267805651450434, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 10.753946960547603, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 10.279264504672039, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 10.688507865341325, 0.0001 * arcsecond)

def test_sun_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['sun']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (1.0160878650466754, 1.0118605934887042, 0.9833276788862821, 1.0107820040799866), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (8.03008088792976, 3.776110727862678, 18.752544254682526, 10.268162490439073), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (20.496475643233936, 19.907832379364574, -23.033309607967187, 10.751933902906119), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (8.02969030304998, 3.775721385487214, 18.752126228091367, 10.267805651450434), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (20.497605463260726, 19.906601181542, -23.03376015263556, 10.753946960547603), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (8.000108116572395, 3.7666292045824337, 18.75183797477899, 10.279264504672039), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (20.58493093599605, 19.879173772309745, -23.032488638722818, 10.688507865341325), 0.0001 * arcsecond)

def test_moon_geocentric_date0(de405, ts):
    t = ts.tt_jd(2440423.345833333)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['moon']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.0026034424248854585, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.472463241145173, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -4.546618838170065, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.472340287066462, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -4.545964408923231, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.446262111681095, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -4.378227942512158, 0.0001 * arcsecond)

def test_moon_geocentric_date1(de405, ts):
    t = ts.tt_jd(2448031.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['moon']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.0024815092296598847, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 23.676443817409496, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 1.8587554901327035, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.676289920709518, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 1.857413875990142, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.66827809687387, 0.0002 * ra_arcsecond)
    compare(dec.degrees, 1.8051891857266409, 0.0001 * arcsecond)

def test_moon_geocentric_date2(de405, ts):
    t = ts.tt_jd(2451545.0)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['moon']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.002690202988513297, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 14.830020573942235, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -10.900635500943373, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 14.829807890359675, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -10.90012775884129, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 14.829573271760747, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -10.897905576904787, 0.0001 * arcsecond)

def test_moon_geocentric_date3(de405, ts):
    t = ts.tt_jd(2456164.5)
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['moon']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, 0.0024739078649309238, 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 16.39102815233177, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -20.93676001523414, 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.39106196861365, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -20.936774891979848, 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.40383113143219, 0.0002 * ra_arcsecond)
    compare(dec.degrees, -20.96508913558473, 0.0001 * arcsecond)

def test_moon_geocentric_date4(de405, ts):
    t = ts.tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    reduce_precision(t)

    e = de405['earth'].at(t)
    p = de405['moon']

    distance = length_of((e - p.at(t)).position.au)
    compare(distance * OLD_AU, (0.0026034424248854585, 0.0024815092296598847, 0.002690202988513297, 0.0024739078649309238), 0.014 * meter)

    astrometric = e.observe(p)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.472463241145173, 23.676443817409496, 14.830020573942235, 16.39102815233177), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-4.546618838170065, 1.8587554901327035, -10.900635500943373, -20.93676001523414), 0.0001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.472340287066462, 23.676289920709518, 14.829807890359675, 16.39106196861365), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-4.545964408923231, 1.857413875990142, -10.90012775884129, -20.936774891979848), 0.0001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.446262111681095, 23.66827809687387, 14.829573271760747, 16.40383113143219), 0.0002 * ra_arcsecond)
    compare(dec.degrees, (-4.378227942512158, 1.8051891857266409, -10.897905576904787, -20.96508913558473), 0.0001 * arcsecond)

def test_polaris_geocentric_date0(earth):
    e = earth.at(load.timescale().tt_jd(2440423.345833333))
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.5283697499529345, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.2642084845529, 0.00001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.52280149297809, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.25882879505869, 0.00001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.0385816433557173, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.11999387030946, 0.00001 * arcsecond)

def test_polaris_geocentric_date1(earth):
    e = earth.at(load.timescale().tt_jd(2448031.5))
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.5296910275944064, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26413894692217, 0.00001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.503356852811078, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26201007627152, 0.00001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.3329211805288432, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.22082922133737, 0.00001 * arcsecond)

def test_polaris_geocentric_date2(earth):
    e = earth.at(load.timescale().tt_jd(2451545.0))
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.5302921882000127, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26411027119273, 0.00001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.544633215462727, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26917874902797, 0.00001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.5459982729094564, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26697328449004, 0.00001 * arcsecond)

def test_polaris_geocentric_date3(earth):
    e = earth.at(load.timescale().tt_jd(2456164.5))
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.531117065610149, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26406906493733, 0.00001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.541609533735535, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.25923373182651, 0.00001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.8064741334456413, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.3136939266471, 0.00001 * arcsecond)

def test_polaris_geocentric_date4(earth):
    e = earth.at(load.timescale().tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5]))
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (2.5283697499529345, 2.5296910275944064, 2.5302921882000127, 2.531117065610149), 0.00001 * ra_arcsecond)
    compare(dec.degrees, (89.2642084845529, 89.26413894692217, 89.26411027119273, 89.26406906493733), 0.00001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (2.52280149297809, 2.503356852811078, 2.544633215462727, 2.541609533735535), 0.00001 * ra_arcsecond)
    compare(dec.degrees, (89.25882879505869, 89.26201007627152, 89.26917874902797, 89.25923373182651), 0.00001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (2.0385816433557173, 2.3329211805288432, 2.5459982729094564, 2.8064741334456413), 0.00001 * ra_arcsecond)
    compare(dec.degrees, (89.11999387030946, 89.22082922133737, 89.26697328449004, 89.3136939266471), 0.00001 * arcsecond)

def test_mercury_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mercury']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.9049140222444105, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 22.33276016366845, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.874898511438327, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 22.415294637224765, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 46.3212267566032, 0.0005 * arcsecond)
    compare(az.degrees, 262.18590521567705, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 46.33688339908365, 0.0005 * arcsecond)
    compare(az.degrees, 262.18590521567705, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 46.33704240110901, 0.0005 * arcsecond)
    compare(az.degrees, 262.18590521567705, 0.0005 * arcsecond)

def test_mercury_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mercury']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.469959592064856, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 11.24594905426479, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.461508188066483, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 11.205182598299666, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -17.340667089884377, 0.0005 * arcsecond)
    compare(az.degrees, 300.9176579181716, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -17.340667089884377, 0.0005 * arcsecond)
    compare(az.degrees, 300.9176579181716, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -17.340667089884377, 0.0005 * arcsecond)
    compare(az.degrees, 300.9176579181716, 0.0005 * arcsecond)

def test_mercury_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mercury']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.138603904058247, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -24.421550562485436, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.138318996641566, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -24.420003066967503, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -0.12765060376706572, 0.0005 * arcsecond)
    compare(az.degrees, 121.97764361867154, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 0.36890915770104016, 0.0005 * arcsecond)
    compare(az.degrees, 121.97764361867154, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 0.3731892291678349, 0.0005 * arcsecond)
    compare(az.degrees, 121.97764361867154, 0.0005 * arcsecond)

def test_mercury_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mercury']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 9.29546814256182, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 16.68590812465023, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 9.307459135231527, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 16.630243128506475, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -9.116616855755964, 0.0005 * arcsecond)
    compare(az.degrees, 300.1420264373104, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -9.116616855755964, 0.0005 * arcsecond)
    compare(az.degrees, 300.1420264373104, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -9.116616855755964, 0.0005 * arcsecond)
    compare(az.degrees, 300.1420264373104, 0.0005 * arcsecond)

def test_mercury_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mercury']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (7.9049140222444105, 2.469959592064856, 18.138603904058247, 9.29546814256182), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (22.33276016366845, 11.24594905426479, -24.421550562485436, 16.68590812465023), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (7.874898511438327, 2.461508188066483, 18.138318996641566, 9.307459135231527), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (22.415294637224765, 11.205182598299666, -24.420003066967503, 16.630243128506475), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (46.3212267566032, -17.340667089884377, -0.12765060376706572, -9.116616855755964), 0.0005 * arcsecond)
    compare(az.degrees, (262.18590521567705, 300.9176579181716, 121.97764361867154, 300.1420264373104), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (46.33688339908365, -17.340667089884377, 0.36890915770104016, -9.116616855755964), 0.0005 * arcsecond)
    compare(az.degrees, (262.18590521567705, 300.9176579181716, 121.97764361867154, 300.1420264373104), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (46.33704240110901, -17.340667089884377, 0.3731892291678349, -9.116616855755964), 0.0005 * arcsecond)
    compare(az.degrees, (262.18590521567705, 300.9176579181716, 121.97764361867154, 300.1420264373104), 0.0005 * arcsecond)

def test_venus_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['venus']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.9665155792599744, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 20.20866872703497, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.936546062416392, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 20.165161469755127, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 11.152374062990575, 0.0005 * arcsecond)
    compare(az.degrees, 287.0030740239532, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 11.23199275246975, 0.0005 * arcsecond)
    compare(az.degrees, 287.0030740239532, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 11.232796262162083, 0.0005 * arcsecond)
    compare(az.degrees, 287.0030740239532, 0.0005 * arcsecond)

def test_venus_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['venus']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.1614662937271143, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 5.325222585955545, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.1533422187037876, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 5.275819541572404, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -34.134914076462266, 0.0005 * arcsecond)
    compare(az.degrees, 313.64872862118426, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -34.134914076462266, 0.0005 * arcsecond)
    compare(az.degrees, 313.64872862118426, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -34.134914076462266, 0.0005 * arcsecond)
    compare(az.degrees, 313.64872862118426, 0.0005 * arcsecond)

def test_venus_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['venus']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.99311221167692, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -18.45256680288619, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.99286396137589, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -18.450404301558034, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 23.228910604670816, 0.0005 * arcsecond)
    compare(az.degrees, 142.1161398141626, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 23.266773672986005, 0.0005 * arcsecond)
    compare(az.degrees, 142.1161398141626, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 23.267157712313676, 0.0005 * arcsecond)
    compare(az.degrees, 142.1161398141626, 0.0005 * arcsecond)

def test_venus_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['venus']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.175218975921811, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 19.87224931182421, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.187940160922054, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 19.849149573371733, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -24.359995410915445, 0.0005 * arcsecond)
    compare(az.degrees, 327.640588969984, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -24.359995410915445, 0.0005 * arcsecond)
    compare(az.degrees, 327.640588969984, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -24.359995410915445, 0.0005 * arcsecond)
    compare(az.degrees, 327.640588969984, 0.0005 * arcsecond)

def test_venus_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['venus']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (4.9665155792599744, 1.1614662937271143, 15.99311221167692, 7.175218975921811), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (20.20866872703497, 5.325222585955545, -18.45256680288619, 19.87224931182421), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (4.936546062416392, 1.1533422187037876, 15.99286396137589, 7.187940160922054), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (20.165161469755127, 5.275819541572404, -18.450404301558034, 19.849149573371733), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (11.152374062990575, -34.134914076462266, 23.228910604670816, -24.359995410915445), 0.0005 * arcsecond)
    compare(az.degrees, (287.0030740239532, 313.64872862118426, 142.1161398141626, 327.640588969984), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (11.23199275246975, -34.134914076462266, 23.266773672986005, -24.359995410915445), 0.0005 * arcsecond)
    compare(az.degrees, (287.0030740239532, 313.64872862118426, 142.1161398141626, 327.640588969984), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (11.232796262162083, -34.134914076462266, 23.267157712313676, -24.359995410915445), 0.0005 * arcsecond)
    compare(az.degrees, (287.0030740239532, 313.64872862118426, 142.1161398141626, 327.640588969984), 0.0005 * arcsecond)

def test_mars_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mars']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.030112454663165, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -24.130883187697044, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.999737237126766, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -24.048966502229923, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -3.540294697028628, 0.0005 * arcsecond)
    compare(az.degrees, 118.34877634707522, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -3.540294697028628, 0.0005 * arcsecond)
    compare(az.degrees, 118.34877634707522, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -3.540294697028628, 0.0005 * arcsecond)
    compare(az.degrees, 118.34877634707522, 0.0005 * arcsecond)

def test_mars_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mars']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.54486790147113, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -4.883946644223003, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.53682348628842, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -4.936042744435578, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -54.1089628741949, 0.0005 * arcsecond)
    compare(az.degrees, 338.0117138951488, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -54.1089628741949, 0.0005 * arcsecond)
    compare(az.degrees, 338.0117138951488, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -54.1089628741949, 0.0005 * arcsecond)
    compare(az.degrees, 338.0117138951488, 0.0005 * arcsecond)

def test_mars_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mars']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.034740913364253, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -13.182784253332377, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.03447079524992, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -13.183338672731741, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -36.90573266459917, 0.0005 * arcsecond)
    compare(az.degrees, 76.12368450672822, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -36.90573266459917, 0.0005 * arcsecond)
    compare(az.degrees, 76.12368450672822, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -36.90573266459917, 0.0005 * arcsecond)
    compare(az.degrees, 76.12368450672822, 0.0005 * arcsecond)

def test_mars_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mars']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.8940809044733, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -12.122804110106655, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.905664739133574, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -12.185661905051244, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 22.094794272017666, 0.0005 * arcsecond)
    compare(az.degrees, 231.6381663847761, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 22.134776069489533, 0.0005 * arcsecond)
    compare(az.degrees, 231.6381663847761, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 22.135181528743814, 0.0005 * arcsecond)
    compare(az.degrees, 231.6381663847761, 0.0005 * arcsecond)

def test_mars_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['mars']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (16.030112454663165, 23.54486790147113, 22.034740913364253, 13.8940809044733), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (-24.130883187697044, -4.883946644223003, -13.182784253332377, -12.122804110106655), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.999737237126766, 23.53682348628842, 22.03447079524992, 13.905664739133574), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (-24.048966502229923, -4.936042744435578, -13.183338672731741, -12.185661905051244), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (-3.540294697028628, -54.1089628741949, -36.90573266459917, 22.094794272017666), 0.0005 * arcsecond)
    compare(az.degrees, (118.34877634707522, 338.0117138951488, 76.12368450672822, 231.6381663847761), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (-3.540294697028628, -54.1089628741949, -36.90573266459917, 22.134776069489533), 0.0005 * arcsecond)
    compare(az.degrees, (118.34877634707522, 338.0117138951488, 76.12368450672822, 231.6381663847761), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (-3.540294697028628, -54.1089628741949, -36.90573266459917, 22.135181528743814), 0.0005 * arcsecond)
    compare(az.degrees, (118.34877634707522, 338.0117138951488, 76.12368450672822, 231.6381663847761), 0.0005 * arcsecond)

def test_jupiter_barycenter_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['jupiter barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.103946503374884, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 0.6522085918269475, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.077992233588102, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 0.821355893113747, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 49.40651603144681, 0.0005 * arcsecond)
    compare(az.degrees, 156.07088561561997, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 49.42056980196601, 0.0005 * arcsecond)
    compare(az.degrees, 156.07088561561997, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 49.420712533159694, 0.0005 * arcsecond)
    compare(az.degrees, 156.07088561561997, 0.0005 * arcsecond)

def test_jupiter_barycenter_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['jupiter barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 6.764836821339949, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 23.17058790055951, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 6.755365668515656, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 23.18253602996423, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 38.00505126690997, 0.0005 * arcsecond)
    compare(az.degrees, 270.63795554820535, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 38.02600464378366, 0.0005 * arcsecond)
    compare(az.degrees, 270.63795554820535, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 38.02621739324931, 0.0005 * arcsecond)
    compare(az.degrees, 270.63795554820535, 0.0005 * arcsecond)

def test_jupiter_barycenter_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['jupiter barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.5914118935512866, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 8.595923929888196, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.5911839414385696, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 8.593862752942394, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -42.482560972481394, 0.0005 * arcsecond)
    compare(az.degrees, 359.3596746827537, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -42.482560972481394, 0.0005 * arcsecond)
    compare(az.degrees, 359.3596746827537, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -42.482560972481394, 0.0005 * arcsecond)
    compare(az.degrees, 359.3596746827537, 0.0005 * arcsecond)

def test_jupiter_barycenter_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['jupiter barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.82276173655752, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 21.649526689253502, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.835667333191383, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 21.670171438742255, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -29.289013841967986, 0.0005 * arcsecond)
    compare(az.degrees, 4.327425566855523, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -29.289013841967986, 0.0005 * arcsecond)
    compare(az.degrees, 4.327425566855523, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -29.289013841967986, 0.0005 * arcsecond)
    compare(az.degrees, 4.327425566855523, 0.0005 * arcsecond)

def test_jupiter_barycenter_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['jupiter barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.103946503374884, 6.764836821339949, 1.5914118935512866, 4.82276173655752), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (0.6522085918269475, 23.17058790055951, 8.595923929888196, 21.649526689253502), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.077992233588102, 6.755365668515656, 1.5911839414385696, 4.835667333191383), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (0.821355893113747, 23.18253602996423, 8.593862752942394, 21.670171438742255), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (49.40651603144681, 38.00505126690997, -42.482560972481394, -29.289013841967986), 0.0005 * arcsecond)
    compare(az.degrees, (156.07088561561997, 270.63795554820535, 359.3596746827537, 4.327425566855523), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (49.42056980196601, 38.02600464378366, -42.482560972481394, -29.289013841967986), 0.0005 * arcsecond)
    compare(az.degrees, (156.07088561561997, 270.63795554820535, 359.3596746827537, 4.327425566855523), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (49.420712533159694, 38.02621739324931, -42.482560972481394, -29.289013841967986), 0.0005 * arcsecond)
    compare(az.degrees, (156.07088561561997, 270.63795554820535, 359.3596746827537, 4.327425566855523), 0.0005 * arcsecond)

def test_saturn_barycenter_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['saturn barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.4626938858905594, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 12.045561201575383, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.4352742791152338, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 11.911391441362444, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -20.662686940324093, 0.0005 * arcsecond)
    compare(az.degrees, 306.01978569992787, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -20.662686940324093, 0.0005 * arcsecond)
    compare(az.degrees, 306.01978569992787, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -20.662686940324093, 0.0005 * arcsecond)
    compare(az.degrees, 306.01978569992787, 0.0005 * arcsecond)

def test_saturn_barycenter_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['saturn barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.814469727768646, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -20.932928080758664, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.805283998285297, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -20.958246345579155, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -48.93337647838982, 0.0005 * arcsecond)
    compare(az.degrees, 76.8837444919445, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -48.93337647838982, 0.0005 * arcsecond)
    compare(az.degrees, 76.8837444919445, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -48.93337647838982, 0.0005 * arcsecond)
    compare(az.degrees, 76.8837444919445, 0.0005 * arcsecond)

def test_saturn_barycenter_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['saturn barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.5845847757319116, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 12.616768688416162, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.584352575888522, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 12.614560194137907, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -36.501918751911674, 0.0005 * arcsecond)
    compare(az.degrees, 341.22347230453323, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -36.501918751911674, 0.0005 * arcsecond)
    compare(az.degrees, 341.22347230453323, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -36.501918751911674, 0.0005 * arcsecond)
    compare(az.degrees, 341.22347230453323, 0.0005 * arcsecond)

def test_saturn_barycenter_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['saturn barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.628268137367913, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -7.658197329820583, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.639621846921335, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -7.723370683249701, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 22.96675851611188, 0.0005 * arcsecond)
    compare(az.degrees, 238.00627672875672, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 23.005094362956072, 0.0005 * arcsecond)
    compare(az.degrees, 238.00627672875672, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 23.005483182929098, 0.0005 * arcsecond)
    compare(az.degrees, 238.00627672875672, 0.0005 * arcsecond)

def test_saturn_barycenter_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['saturn barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (2.4626938858905594, 19.814469727768646, 2.5845847757319116, 13.628268137367913), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (12.045561201575383, -20.932928080758664, 12.616768688416162, -7.658197329820583), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (2.4352742791152338, 19.805283998285297, 2.584352575888522, 13.639621846921335), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (11.911391441362444, -20.958246345579155, 12.614560194137907, -7.723370683249701), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (-20.662686940324093, -48.93337647838982, -36.501918751911674, 22.96675851611188), 0.0005 * arcsecond)
    compare(az.degrees, (306.01978569992787, 76.8837444919445, 341.22347230453323, 238.00627672875672), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (-20.662686940324093, -48.93337647838982, -36.501918751911674, 23.005094362956072), 0.0005 * arcsecond)
    compare(az.degrees, (306.01978569992787, 76.8837444919445, 341.22347230453323, 238.00627672875672), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (-20.662686940324093, -48.93337647838982, -36.501918751911674, 23.005483182929098), 0.0005 * arcsecond)
    compare(az.degrees, (306.01978569992787, 76.8837444919445, 341.22347230453323, 238.00627672875672), 0.0005 * arcsecond)

def test_uranus_barycenter_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['uranus barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.087016642067397, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 0.20824442104711183, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.061058763070791, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 0.37741883683460087, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 49.06396822144731, 0.0005 * arcsecond)
    compare(az.degrees, 156.65256040205296, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 49.078192535060566, 0.0005 * arcsecond)
    compare(az.degrees, 156.65256040205296, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 49.07833699756142, 0.0005 * arcsecond)
    compare(az.degrees, 156.65256040205296, 0.0005 * arcsecond)

def test_uranus_barycenter_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['uranus barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.668863148648313, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -23.43704804377175, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.65936510933368, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -23.447712978993913, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -37.0259637798912, 0.0005 * arcsecond)
    compare(az.degrees, 91.80748703145906, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -37.0259637798912, 0.0005 * arcsecond)
    compare(az.degrees, 91.80748703145906, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -37.0259637798912, 0.0005 * arcsecond)
    compare(az.degrees, 91.80748703145906, 0.0005 * arcsecond)

def test_uranus_barycenter_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['uranus barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 21.16527335872666, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -17.020308119118386, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 21.164991487815, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -17.020361566142082, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -29.175475562665554, 0.0005 * arcsecond)
    compare(az.degrees, 88.85671230431439, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -29.175475562665554, 0.0005 * arcsecond)
    compare(az.degrees, 88.85671230431439, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -29.175475562665554, 0.0005 * arcsecond)
    compare(az.degrees, 88.85671230431439, 0.0005 * arcsecond)

def test_uranus_barycenter_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['uranus barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 0.48945083888242796, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 2.358286196725548, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 0.5005545778924997, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 2.4296958868419787, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -14.5260443119261, 0.0005 * arcsecond)
    compare(az.degrees, 74.60219420538265, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -14.5260443119261, 0.0005 * arcsecond)
    compare(az.degrees, 74.60219420538265, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -14.5260443119261, 0.0005 * arcsecond)
    compare(az.degrees, 74.60219420538265, 0.0005 * arcsecond)

def test_uranus_barycenter_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['uranus barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.087016642067397, 18.668863148648313, 21.16527335872666, 0.48945083888242796), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (0.20824442104711183, -23.43704804377175, -17.020308119118386, 2.358286196725548), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.061058763070791, 18.65936510933368, 21.164991487815, 0.5005545778924997), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (0.37741883683460087, -23.447712978993913, -17.020361566142082, 2.4296958868419787), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (49.06396822144731, -37.0259637798912, -29.175475562665554, -14.5260443119261), 0.0005 * arcsecond)
    compare(az.degrees, (156.65256040205296, 91.80748703145906, 88.85671230431439, 74.60219420538265), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (49.078192535060566, -37.0259637798912, -29.175475562665554, -14.5260443119261), 0.0005 * arcsecond)
    compare(az.degrees, (156.65256040205296, 91.80748703145906, 88.85671230431439, 74.60219420538265), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (49.07833699756142, -37.0259637798912, -29.175475562665554, -14.5260443119261), 0.0005 * arcsecond)
    compare(az.degrees, (156.65256040205296, 91.80748703145906, 88.85671230431439, 74.60219420538265), 0.0005 * arcsecond)

def test_neptune_barycenter_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['neptune barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.637396931781986, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -17.680489951171502, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.608492665044128, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -17.583829722494027, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 4.86937782636538, 0.0005 * arcsecond)
    compare(az.degrees, 117.29043762875409, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 5.031511017145419, 0.0005 * arcsecond)
    compare(az.degrees, 117.29043762875409, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 5.033116634143141, 0.0005 * arcsecond)
    compare(az.degrees, 117.29043762875409, 0.0005 * arcsecond)

def test_neptune_barycenter_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['neptune barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.036514568239326, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -21.792523874854822, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.027165016434417, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -21.808061138689617, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -40.43318694811052, 0.0005 * arcsecond)
    compare(az.degrees, 86.51833613444356, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -40.43318694811052, 0.0005 * arcsecond)
    compare(az.degrees, 86.51833613444356, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -40.43318694811052, 0.0005 * arcsecond)
    compare(az.degrees, 86.51833613444356, 0.0005 * arcsecond)

def test_neptune_barycenter_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['neptune barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 20.362478654099593, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -19.213665913911328, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 20.36219137258442, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -19.21325376377245, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -21.102154672787563, 0.0005 * arcsecond)
    compare(az.degrees, 98.14962081515444, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -21.102154672787563, 0.0005 * arcsecond)
    compare(az.degrees, 98.14962081515444, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -21.102154672787563, 0.0005 * arcsecond)
    compare(az.degrees, 98.14962081515444, 0.0005 * arcsecond)

def test_neptune_barycenter_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['neptune barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.252831344843074, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -11.502690543226894, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.26432121506238, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -11.437371208596403, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 2.41678290499992, 0.0005 * arcsecond)
    compare(az.degrees, 106.8092597257607, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 2.6713913487620147, 0.0005 * arcsecond)
    compare(az.degrees, 106.8092597257607, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 2.6738334093305696, 0.0005 * arcsecond)
    compare(az.degrees, 106.8092597257607, 0.0005 * arcsecond)

def test_neptune_barycenter_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['neptune barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (15.637396931781986, 19.036514568239326, 20.362478654099593, 22.252831344843074), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (-17.680489951171502, -21.792523874854822, -19.213665913911328, -11.502690543226894), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.608492665044128, 19.027165016434417, 20.36219137258442, 22.26432121506238), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (-17.583829722494027, -21.808061138689617, -19.21325376377245, -11.437371208596403), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (4.86937782636538, -40.43318694811052, -21.102154672787563, 2.41678290499992), 0.0005 * arcsecond)
    compare(az.degrees, (117.29043762875409, 86.51833613444356, 98.14962081515444, 106.8092597257607), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (5.031511017145419, -40.43318694811052, -21.102154672787563, 2.6713913487620147), 0.0005 * arcsecond)
    compare(az.degrees, (117.29043762875409, 86.51833613444356, 98.14962081515444, 106.8092597257607), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (5.033116634143141, -40.43318694811052, -21.102154672787563, 2.6738334093305696), 0.0005 * arcsecond)
    compare(az.degrees, (117.29043762875409, 86.51833613444356, 98.14962081515444, 106.8092597257607), 0.0005 * arcsecond)

def test_pluto_barycenter_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['pluto barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.015146948702718, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 16.622956629676764, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 11.989238323883423, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 16.792209116103148, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 64.72856074651983, 0.0005 * arcsecond)
    compare(az.degrees, 147.2138070056058, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 64.73630449169308, 0.0005 * arcsecond)
    compare(az.degrees, 147.2138070056058, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 64.73638314930092, 0.0005 * arcsecond)
    compare(az.degrees, 147.2138070056058, 0.0005 * arcsecond)

def test_pluto_barycenter_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['pluto barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.216666873470118, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -1.335915234746897, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.208587498665665, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -1.3022917220648205, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 16.233734452123414, 0.0005 * arcsecond)
    compare(az.degrees, 105.3994365631196, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 16.28889280191291, 0.0005 * arcsecond)
    compare(az.degrees, 105.3994365631196, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 16.289451329649054, 0.0005 * arcsecond)
    compare(az.degrees, 105.3994365631196, 0.0005 * arcsecond)

def test_pluto_barycenter_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['pluto barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.761532920101487, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -11.396347593297179, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.76128368305737, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -11.39433478419375, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 22.700996363632996, 0.0005 * arcsecond)
    compare(az.degrees, 127.81134408260581, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 22.739821647292274, 0.0005 * arcsecond)
    compare(az.degrees, 127.81134408260581, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 22.74021541578692, 0.0005 * arcsecond)
    compare(az.degrees, 127.81134408260581, 0.0005 * arcsecond)

def test_pluto_barycenter_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['pluto barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.488579709427018, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -19.551785355075808, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.501344365322606, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -19.541283736216652, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 28.33982206878914, 0.0005 * arcsecond)
    compare(az.degrees, 157.51785266272373, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 28.370071242061236, 0.0005 * arcsecond)
    compare(az.degrees, 157.51785266272373, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 28.370378222043662, 0.0005 * arcsecond)
    compare(az.degrees, 157.51785266272373, 0.0005 * arcsecond)

def test_pluto_barycenter_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['pluto barycenter']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.015146948702718, 15.216666873470118, 16.761532920101487, 18.488579709427018), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (16.622956629676764, -1.335915234746897, -11.396347593297179, -19.551785355075808), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (11.989238323883423, 15.208587498665665, 16.76128368305737, 18.501344365322606), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (16.792209116103148, -1.3022917220648205, -11.39433478419375, -19.541283736216652), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (64.72856074651983, 16.233734452123414, 22.700996363632996, 28.33982206878914), 0.0005 * arcsecond)
    compare(az.degrees, (147.2138070056058, 105.3994365631196, 127.81134408260581, 157.51785266272373), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (64.73630449169308, 16.28889280191291, 22.739821647292274, 28.370071242061236), 0.0005 * arcsecond)
    compare(az.degrees, (147.2138070056058, 105.3994365631196, 127.81134408260581, 157.51785266272373), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (64.73638314930092, 16.289451329649054, 22.74021541578692, 28.370378222043662), 0.0005 * arcsecond)
    compare(az.degrees, (147.2138070056058, 105.3994365631196, 127.81134408260581, 157.51785266272373), 0.0005 * arcsecond)

def test_sun_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['sun']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 8.02959789881544, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 20.496678572125123, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 8.000015838288707, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 20.584000539289498, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 46.72403357148823, 0.0005 * arcsecond)
    compare(az.degrees, 258.5550717845957, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 46.73947196634687, 0.0005 * arcsecond)
    compare(az.degrees, 258.5550717845957, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 46.73962875307724, 0.0005 * arcsecond)
    compare(az.degrees, 258.5550717845957, 0.0005 * arcsecond)

def test_sun_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['sun']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 3.7755906381611175, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 19.90505409109931, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 3.7664985705990794, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 19.87762515818775, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 2.2209469369832533, 0.0005 * arcsecond)
    compare(az.degrees, 293.95636637272145, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 2.4868409787793837, 0.0005 * arcsecond)
    compare(az.degrees, 293.95636637272145, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 2.489379891081029, 0.0005 * arcsecond)
    compare(az.degrees, 293.95636637272145, 0.0005 * arcsecond)

def test_sun_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['sun']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.752264357691004, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -23.03532101826747, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.751976099155204, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -23.03404957045815, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -5.486505415022805, 0.0005 * arcsecond)
    compare(az.degrees, 115.32008451470392, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -5.486505415022805, 0.0005 * arcsecond)
    compare(az.degrees, 115.32008451470392, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -5.486505415022805, 0.0005 * arcsecond)
    compare(az.degrees, 115.32008451470392, 0.0005 * arcsecond)

def test_sun_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['sun']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 10.267679924967121, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 10.752399537108259, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 10.279138748598198, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 10.686961444410377, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -2.738407691502772, 0.0005 * arcsecond)
    compare(az.degrees, 286.09632001391725, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -2.738407691502772, 0.0005 * arcsecond)
    compare(az.degrees, 286.09632001391725, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -2.738407691502772, 0.0005 * arcsecond)
    compare(az.degrees, 286.09632001391725, 0.0005 * arcsecond)

def test_sun_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['sun']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (8.02959789881544, 3.7755906381611175, 18.752264357691004, 10.267679924967121), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (20.496678572125123, 19.90505409109931, -23.03532101826747, 10.752399537108259), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (8.000015838288707, 3.7664985705990794, 18.751976099155204, 10.279138748598198), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (20.584000539289498, 19.87762515818775, -23.03404957045815, 10.686961444410377), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (46.72403357148823, 2.2209469369832533, -5.486505415022805, -2.738407691502772), 0.0005 * arcsecond)
    compare(az.degrees, (258.5550717845957, 293.95636637272145, 115.32008451470392, 286.09632001391725), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (46.73947196634687, 2.4868409787793837, -5.486505415022805, -2.738407691502772), 0.0005 * arcsecond)
    compare(az.degrees, (258.5550717845957, 293.95636637272145, 115.32008451470392, 286.09632001391725), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (46.73962875307724, 2.489379891081029, -5.486505415022805, -2.738407691502772), 0.0005 * arcsecond)
    compare(az.degrees, (258.5550717845957, 293.95636637272145, 115.32008451470392, 286.09632001391725), 0.0005 * arcsecond)

def test_moon_topocentric_date0(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2440423.345833333)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['moon']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.489955349304845, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -5.189705732227236, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.463855411284248, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -5.022075882872161, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 41.92040135025528, 0.0005 * arcsecond)
    compare(az.degrees, 151.19707488767745, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 41.938650930940234, 0.0005 * arcsecond)
    compare(az.degrees, 151.19707488767745, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 41.938836248377605, 0.0005 * arcsecond)
    compare(az.degrees, 151.19707488767745, 0.0005 * arcsecond)

def test_moon_topocentric_date1(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2448031.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['moon']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.663473338211578, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 1.227161288913488, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.655459675858083, 0.0005 * ra_arcsecond)
    compare(dec.degrees, 1.1749464194383863, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -47.74510120858602, 0.0005 * arcsecond)
    compare(az.degrees, 338.13295291812307, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -47.74510120858602, 0.0005 * arcsecond)
    compare(az.degrees, 338.13295291812307, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -47.74510120858602, 0.0005 * arcsecond)
    compare(az.degrees, 338.13295291812307, 0.0005 * arcsecond)

def test_moon_topocentric_date2(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2451545.0)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['moon']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 14.845679251156893, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -11.590214641232205, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 14.845444624832663, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -11.58799188846256, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 36.381265580736255, 0.0005 * arcsecond)
    compare(az.degrees, 156.2971102404744, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 36.40348032108563, 0.0005 * arcsecond)
    compare(az.degrees, 156.2971102404744, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 36.403705864717445, 0.0005 * arcsecond)
    compare(az.degrees, 156.2971102404744, 0.0005 * arcsecond)

def test_moon_topocentric_date3(de405):
    t = load.timescale(delta_t=0.0).tt_jd(2456164.5)
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['moon']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.380804513901573, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -21.79048462924397, 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.393647715389825, 0.0005 * ra_arcsecond)
    compare(dec.degrees, -21.81897641768761, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 28.439387966372543, 0.0005 * arcsecond)
    compare(az.degrees, 191.29497427201525, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 28.46951344291743, 0.0005 * arcsecond)
    compare(az.degrees, 191.29497427201525, 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 28.46981916998486, 0.0005 * arcsecond)
    compare(az.degrees, 191.29497427201525, 0.0005 * arcsecond)

def test_moon_topocentric_date4(de405):
    t = load.timescale(delta_t=0.0).tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    earth = de405['earth']
    usno = earth + Topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno.at(t).observe(de405['moon']).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.489955349304845, 23.663473338211578, 14.845679251156893, 16.380804513901573), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (-5.189705732227236, 1.227161288913488, -11.590214641232205, -21.79048462924397), 0.0005 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.463855411284248, 23.655459675858083, 14.845444624832663, 16.393647715389825), 0.0005 * ra_arcsecond)
    compare(dec.degrees, (-5.022075882872161, 1.1749464194383863, -11.58799188846256, -21.81897641768761), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (41.92040135025528, -47.74510120858602, 36.381265580736255, 28.439387966372543), 0.0005 * arcsecond)
    compare(az.degrees, (151.19707488767745, 338.13295291812307, 156.2971102404744, 191.29497427201525), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (41.938650930940234, -47.74510120858602, 36.40348032108563, 28.46951344291743), 0.0005 * arcsecond)
    compare(az.degrees, (151.19707488767745, 338.13295291812307, 156.2971102404744, 191.29497427201525), 0.0005 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (41.938836248377605, -47.74510120858602, 36.403705864717445, 28.46981916998486), 0.0005 * arcsecond)
    compare(az.degrees, (151.19707488767745, 338.13295291812307, 156.2971102404744, 191.29497427201525), 0.0005 * arcsecond)

def test_hipparcos_conversion0(earth):
    line = b'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    df = hipparcos.load_dataframe(BytesIO(line))
    star = starlib.Star.from_dataframe(df.iloc[0])
    ra, dec, distance = earth.at(load.timescale().tt_jd(2440423.345833333)).observe(star).radec()
    compare(ra.hours, 2.5283697000528966, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26420852419295, 0.00001 * arcsecond)

def test_hipparcos_conversion1(earth):
    line = b'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    df = hipparcos.load_dataframe(BytesIO(line))
    star = starlib.Star.from_dataframe(df.iloc[0])
    ra, dec, distance = earth.at(load.timescale().tt_jd(2448031.5)).observe(star).radec()
    compare(ra.hours, 2.529691010447949, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26413900274704, 0.00001 * arcsecond)

def test_hipparcos_conversion2(earth):
    line = b'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    df = hipparcos.load_dataframe(BytesIO(line))
    star = starlib.Star.from_dataframe(df.iloc[0])
    ra, dec, distance = earth.at(load.timescale().tt_jd(2451545.0)).observe(star).radec()
    compare(ra.hours, 2.5302921836971946, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26411033462212, 0.00001 * arcsecond)

def test_hipparcos_conversion3(earth):
    line = b'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    df = hipparcos.load_dataframe(BytesIO(line))
    star = starlib.Star.from_dataframe(df.iloc[0])
    ra, dec, distance = earth.at(load.timescale().tt_jd(2456164.5)).observe(star).radec()
    compare(ra.hours, 2.5311170753257395, 0.00001 * ra_arcsecond)
    compare(dec.degrees, 89.26406913848278, 0.00001 * arcsecond)

def test_hipparcos_conversion4(earth):
    line = b'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    df = hipparcos.load_dataframe(BytesIO(line))
    star = starlib.Star.from_dataframe(df.iloc[0])
    ra, dec, distance = earth.at(load.timescale().tt_jd([2440423.345833333, 2448031.5, 2451545.0, 2456164.5])).observe(star).radec()
    compare(ra.hours, (2.5283697000528966, 2.529691010447949, 2.5302921836971946, 2.5311170753257395), 0.00001 * ra_arcsecond)
    compare(dec.degrees, (89.26420852419295, 89.26413900274704, 89.26411033462212, 89.26406913848278), 0.00001 * arcsecond)

