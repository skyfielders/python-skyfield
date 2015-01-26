'Auto-generated accuracy tests vs NOVAS (see build_novas_tests.py).'

import assay
from numpy import abs, array, einsum, max
from skyfield import (earthlib, framelib, nutationlib, positionlib,
                      precessionlib, starlib, timelib)
from skyfield.api import JulianDate
from skyfield.constants import AU_M
from skyfield.data import hipparcos
from skyfield.functions import length_of
from skyfield.jpllib import Ephemeris
from skyfield.positionlib import Apparent

import de405
de405 = Ephemeris(de405)

one_second = 1.0 / 24.0 / 60.0 / 60.0
arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
meter = 1.0 / AU_M

def compare(value, expected_value, epsilon):
    if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
        assert max(abs(value - expected_value)) <= epsilon
    else:
        assert abs(value - expected_value) <= epsilon

def test_calendar_date_0():
    compare(timelib.calendar_date(2440423.345833333), array((1969, 7, 20.345833333209157)), 0.0)

def test_calendar_date_1():
    compare(timelib.calendar_date(2448031.5), array((1990, 5, 19.5)), 0.0)

def test_calendar_date_2():
    compare(timelib.calendar_date(2451545.0), array((2000, 1, 1.0)), 0.0)

def test_calendar_date_3():
    compare(timelib.calendar_date(2456164.5), array((2012, 8, 24.5)), 0.0)

def test_earth_rotation_angle_date0():
    compare(earthlib.earth_rotation_angle(2440423.345833333) * 360.0, 243.32160780274833,
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

def test_earth_tilt_date0():
    compare(nutationlib.earth_tilt(JulianDate(tdb=2440423.345833333)),
            array((23.443240959852666, 23.445702723464045, 0.1592945569695419, 2.604727521416371, 8.862349000962688)), 0.00001 * arcsecond)

def test_earth_tilt_date1():
    compare(nutationlib.earth_tilt(JulianDate(tdb=2448031.5)),
            array((23.440530953006782, 23.442178709915066, 0.7110982205507697, 11.62814814196408, 5.931924869819573)), 0.00001 * arcsecond)

def test_earth_tilt_date2():
    compare(nutationlib.earth_tilt(JulianDate(tdb=2451545.0)),
            array((23.439279444444445, 23.437676833867652, -0.8520167470908031, -13.931996330960066, -5.769398076465292)), 0.00001 * arcsecond)

def test_earth_tilt_date3():
    compare(nutationlib.earth_tilt(JulianDate(tdb=2456164.5)),
            array((23.43763397776759, 23.43645066577372, 0.977087608170225, 15.976729533480201, -4.259923177932991)), 0.00001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date0():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2440423.345833333),
            array(-1.4592438843164885e-09), 0.0000000000000001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date1():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2448031.5),
            array(-9.909270679336256e-09), 0.0000000000000001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date2():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2451545.0),
            array(1.0213300963024648e-08), 0.0000000000000001 * arcsecond)

def test_equation_of_the_equinoxes_complimentary_terms_date3():
    compare(nutationlib.equation_of_the_equinoxes_complimentary_terms(2456164.5),
            array(-1.082315527387237e-08), 0.0000000000000001 * arcsecond)

def test_forward_frame_tie():
    compare(framelib.ICRS_to_J2000.dot((1.1, 1.2, 1.3)), (1.100000019790573, 1.2000001208396125, 1.2999998717098593), 1e-15)

def test_reverse_frame_tie():
    compare(framelib.ICRS_to_J2000.T.dot((1.1, 1.2, 1.3)), (1.0999999802094145, 1.1999998791603803, 1.300000128290131), 1e-15)

def test_fundamental_arguments_date0():
    compare(nutationlib.fundamental_arguments(-0.3044942961441969),
            array((-1.5597846169353031, -2.8619278194907483, -2.7748368269156427, -4.947060102171707, 6.178085194718492)), 0.000000002 * arcsecond)

def test_fundamental_arguments_date1():
    compare(nutationlib.fundamental_arguments(-0.09619438740588637),
            array((-0.8532784044768771, -3.933579124091533, -5.376486844354975, -0.9485312704748627, 5.429677887938805)), 0.000000002 * arcsecond)

def test_fundamental_arguments_date2():
    compare(nutationlib.fundamental_arguments(0.0),
            array((2.355555743493879, 6.24006012692298, 1.6279050815375191, 5.198466588650503, 2.182439196615671)), 0.000000002 * arcsecond)

def test_fundamental_arguments_date3():
    compare(nutationlib.fundamental_arguments(0.12647501711156742),
            array((0.15181719486225662, 4.0231516222224535, 0.10917837795923366, 1.6234303368860354, -2.086983188457769)), 0.000000002 * arcsecond)

def test_iau2000a_date0():
    compare(nutationlib.iau2000a(2440423.345833333),
            array([26047275.21416371, 88623490.00962688]), 0.001)

def test_iau2000a_date1():
    compare(nutationlib.iau2000a(2448031.5),
            array([116281481.41964081, 59319248.69819573]), 0.001)

def test_iau2000a_date2():
    compare(nutationlib.iau2000a(2451545.0),
            array([-139319963.30960065, -57693980.76465292]), 0.001)

def test_iau2000a_date3():
    compare(nutationlib.iau2000a(2456164.5),
            array([159767295.334802, -42599231.77932991]), 0.001)

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

def test_nutation_date0():
    matrix = nutationlib.compute_nutation(JulianDate(tdb=2440423.345833333))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0999795659425045, 1.1999568871469581, 1.3000570847072532),
            result, 1e-14)

def test_nutation_date1():
    matrix = nutationlib.compute_nutation(JulianDate(tdb=2448031.5))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0999087778623433, 1.2000195046911977, 1.300059178938428),
            result, 1e-14)

def test_nutation_date2():
    matrix = nutationlib.compute_nutation(JulianDate(tdb=2451545.0))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.100109290032102, 1.1999681897164485, 1.2999368806421698),
            result, 1e-14)

def test_nutation_date3():
    matrix = nutationlib.compute_nutation(JulianDate(tdb=2456164.5))
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0998746654010052, 1.2001050177909844, 1.3000091025381042),
            result, 1e-14)

def test_precession_date0():
    matrix = precessionlib.compute_precession(2440423.345833333)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.111985657355239, 1.1924703352076302, 1.296727572578774),
            result, 1e-15)

def test_precession_date1():
    matrix = precessionlib.compute_precession(2448031.5)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.103793140541002, 1.1976299348492718, 1.2989700697273823),
            result, 1e-15)

def test_precession_date2():
    matrix = precessionlib.compute_precession(2451545.0)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.1, 1.2, 1.3),
            result, 1e-15)

def test_precession_date3():
    matrix = precessionlib.compute_precession(2456164.5)
    result = einsum('ij...,j...->i...', matrix, [1.1, 1.2, 1.3])
    compare((1.0950034772583117, 1.2031039092689229, 1.301348672836777),
            result, 1e-15)

def test_sidereal_time_on_date0():
    jd = JulianDate(tt=2440423.345833333)
    compare(earthlib.sidereal_time(jd), 16.19543622705723, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date0():
    jd = JulianDate(tt=2440423.345833333 + 99.9 * one_second, delta_t=99.9)
    compare(earthlib.sidereal_time(jd), 16.195436229760517, 1e-13)

def test_sidereal_time_on_date1():
    jd = JulianDate(tt=2448031.5)
    compare(earthlib.sidereal_time(jd), 15.825907460288224, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date1():
    jd = JulianDate(tt=2448031.5 + 99.9 * one_second, delta_t=99.9)
    compare(earthlib.sidereal_time(jd), 15.825907462991848, 1e-13)

def test_sidereal_time_on_date2():
    jd = JulianDate(tt=2451545.0)
    compare(earthlib.sidereal_time(jd), 18.69737482696563, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date2():
    jd = JulianDate(tt=2451545.0 + 99.9 * one_second, delta_t=99.9)
    compare(earthlib.sidereal_time(jd), 18.69737482966941, 1e-13)

def test_sidereal_time_on_date3():
    jd = JulianDate(tt=2456164.5)
    compare(earthlib.sidereal_time(jd), 22.243908497165812, 1e-13)

def test_sidereal_time_with_nonzero_delta_t_on_date3():
    jd = JulianDate(tt=2456164.5 + 99.9 * one_second, delta_t=99.9)
    compare(earthlib.sidereal_time(jd), 22.2439084998698, 1e-13)

def test_star_vector():
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)
    compare(star._position_AU,
            (276301.5236796437, 215517.39549460335, 27281454.187831223),
            1e3 * meter)
    compare(star._velocity_AU_per_d,
            (-0.006595734315371152, 0.015163885823867606, -0.010102577482634966),
            1e-6 * meter)

def test_refraction0():
    r = earthlib.refraction(-5, 10, 1010)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction1():
    r = earthlib.refraction(-5, 10, 1013.25)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction2():
    r = earthlib.refraction(-5, 25, 1010)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction3():
    r = earthlib.refraction(-5, 25, 1013.25)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction4():
    r = earthlib.refraction(-1, 10, 1010)
    compare(r, 0.8296919418249878, 0.001 * arcsecond)

def test_refraction5():
    r = earthlib.refraction(-1, 10, 1013.25)
    compare(r, 0.8323617426278901, 0.001 * arcsecond)

def test_refraction6():
    r = earthlib.refraction(-1, 25, 1010)
    compare(r, 0.7879289246190321, 0.001 * arcsecond)

def test_refraction7():
    r = earthlib.refraction(-1, 25, 1013.25)
    compare(r, 0.7904643394754794, 0.001 * arcsecond)

def test_refraction8():
    r = earthlib.refraction(15, 10, 1010)
    compare(r, 0.06056215494995108, 0.001 * arcsecond)

def test_refraction9():
    r = earthlib.refraction(15, 10, 1013.25)
    compare(r, 0.06075703317132469, 0.001 * arcsecond)

def test_refraction10():
    r = earthlib.refraction(15, 25, 1010)
    compare(r, 0.057513724331664955, 0.001 * arcsecond)

def test_refraction11():
    r = earthlib.refraction(15, 25, 1013.25)
    compare(r, 0.057698793246593584, 0.001 * arcsecond)

def test_refraction12():
    r = earthlib.refraction(89.95, 10, 1010)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction13():
    r = earthlib.refraction(89.95, 10, 1013.25)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction14():
    r = earthlib.refraction(89.95, 25, 1010)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refraction15():
    r = earthlib.refraction(89.95, 25, 1013.25)
    compare(r, 0.0, 0.001 * arcsecond)

def test_refract0():
    alt = earthlib.refract(-90, 10.0, 1010.0)
    compare(alt, -90.0, 0.000000001 * arcsecond)

def test_refract1():
    alt = earthlib.refract(-2, 10.0, 1010.0)
    compare(alt, -2.0, 0.000000001 * arcsecond)

def test_refract2():
    alt = earthlib.refract(-1, 10.0, 1010.0)
    compare(alt, -0.34540033564054795, 0.000000001 * arcsecond)

def test_refract3():
    alt = earthlib.refract(0, 10.0, 1010.0)
    compare(alt, 0.4819388815393779, 0.000000001 * arcsecond)

def test_refract4():
    alt = earthlib.refract(1, 10.0, 1010.0)
    compare(alt, 1.362447444478633, 0.000000001 * arcsecond)

def test_refract5():
    alt = earthlib.refract(3, 10.0, 1010.0)
    compare(alt, 3.227564692764261, 0.000000001 * arcsecond)

def test_refract6():
    alt = earthlib.refract(9, 10.0, 1010.0)
    compare(alt, 9.098059272393698, 0.000000001 * arcsecond)

def test_refract7():
    alt = earthlib.refract(90, 10.0, 1010.0)
    compare(alt, 90.0, 0.000000001 * arcsecond)

def test_ITRF_to_GCRS_conversion_on_date0():
    jd = JulianDate(tt=2440423.345833333, delta_t=39.707)
    position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (0.5701172053657788, -1.5232987806096516, 1.3017400651201705), 1e-13)

def test_ITRF_to_GCRS_conversion_on_date1():
    jd = JulianDate(tt=2448031.5, delta_t=57.1136)
    position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (0.4136264927956479, -1.574108193365244, 1.3004216700893525), 1e-13)

def test_ITRF_to_GCRS_conversion_on_date2():
    jd = JulianDate(tt=2451545.0, delta_t=63.8285)
    position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (1.375700857396341, -0.8702954291925735, 1.3000126987400913), 1e-13)

def test_ITRF_to_GCRS_conversion_on_date3():
    jd = JulianDate(tt=2456164.5, delta_t=66.7846)
    position = positionlib.ITRF_to_GCRS(jd, [1.1, 1.2, 1.3])
    compare(position, (1.5243574049688486, 0.5755748855663745, 1.2980940077752077), 1e-13)

def test_tdb_minus_tt_on_date0():
    result = timelib.tdb_minus_tt(2440423.345833333)
    compare(result, -0.00046798717637515125, 1e-16)

def test_tdb_minus_tt_on_date1():
    result = timelib.tdb_minus_tt(2448031.5)
    compare(result, 0.001158518592634921, 1e-16)

def test_tdb_minus_tt_on_date2():
    result = timelib.tdb_minus_tt(2451545.0)
    compare(result, -9.575743486095212e-05, 1e-16)

def test_tdb_minus_tt_on_date3():
    result = timelib.tdb_minus_tt(2456164.5)
    compare(result, -0.001241030165936061, 1e-16)

def test_mercury_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 1.3278115470600746, 0.5 * meter)

    astrometric = e.observe(de405.mercury)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 7.905384000977572, 0.001 * ra_arcsecond)
    compare(dec.degrees, 22.332364359841474, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.904987228126012, 0.001 * ra_arcsecond)
    compare(dec.degrees, 22.33343308790883, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.874971625095716, 0.001 * ra_arcsecond)
    compare(dec.degrees, 22.41597039204466, 0.001 * arcsecond)

def test_mercury_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 0.6507044512046538, 0.5 * meter)

    astrometric = e.observe(de405.mercury)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.4704717994133576, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.250132844930501, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.470128253572967, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.248550502940757, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.461676722646477, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.20778549324496, 0.001 * arcsecond)

def test_mercury_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 1.4155249674526948, 0.5 * meter)

    astrometric = e.observe(de405.mercury)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.13892977357885, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.420324941080732, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.13851035907211, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.420393338459686, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.138225455402914, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.418845803732097, 0.001 * arcsecond)

def test_mercury_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 1.1264323486728112, 0.5 * meter)

    astrometric = e.observe(de405.mercury)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 9.295934662566733, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.68579742896488, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 9.295575039086721, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.68740973196494, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 9.307566088097712, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.631743449679668, 0.001 * arcsecond)

def test_mercury_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, (1.3278115470600746, 0.6507044512046538, 1.4155249674526948, 1.1264323486728112), 0.5 * meter)

    astrometric = e.observe(de405.mercury)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (7.905384000977572, 2.4704717994133576, 18.13892977357885, 9.295934662566733), 0.001 * ra_arcsecond)
    compare(dec.degrees, (22.332364359841474, 11.250132844930501, -24.420324941080732, 16.68579742896488), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (7.904987228126012, 2.470128253572967, 18.13851035907211, 9.295575039086721), 0.001 * ra_arcsecond)
    compare(dec.degrees, (22.33343308790883, 11.248550502940757, -24.420393338459686, 16.68740973196494), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (7.874971625095716, 2.461676722646477, 18.138225455402914, 9.307566088097712), 0.001 * ra_arcsecond)
    compare(dec.degrees, (22.41597039204466, 11.20778549324496, -24.418845803732097, 16.631743449679668), 0.001 * arcsecond)

def test_venus_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 0.9646045654448725, 0.5 * meter)

    astrometric = e.observe(de405.venus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 4.966946050917652, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.210417323471006, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.966656139420439, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.21014591709747, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.93668626355443, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.166644671858098, 0.001 * arcsecond)

def test_venus_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 1.0711674186789975, 0.5 * meter)

    astrometric = e.observe(de405.venus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 1.161811406279447, 0.001 * ra_arcsecond)
    compare(dec.degrees, 5.32829157368082, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.161541590682067, 0.001 * ra_arcsecond)
    compare(dec.degrees, 5.326768071513866, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.1534174784892792, 0.001 * ra_arcsecond)
    compare(dec.degrees, 5.277365365528824, 0.001 * arcsecond)

def test_venus_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 1.1376890757925104, 0.5 * meter)

    astrometric = e.observe(de405.venus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 15.993350650200568, 0.001 * ra_arcsecond)
    compare(dec.degrees, -18.451653207795232, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.993038357924485, 0.001 * ra_arcsecond)
    compare(dec.degrees, -18.450881488018126, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.99279010971033, 0.001 * ra_arcsecond)
    compare(dec.degrees, -18.448718976425834, 0.001 * arcsecond)

def test_venus_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 0.7824924286112767, 0.5 * meter)

    astrometric = e.observe(de405.venus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 7.175585125577371, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.874130272238094, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.175312328808405, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.874779975491407, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.188033727750362, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.851678563902254, 0.001 * arcsecond)

def test_venus_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, (0.9646045654448725, 1.0711674186789975, 1.1376890757925104, 0.7824924286112767), 0.5 * meter)

    astrometric = e.observe(de405.venus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (4.966946050917652, 1.161811406279447, 15.993350650200568, 7.175585125577371), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.210417323471006, 5.32829157368082, -18.451653207795232, 19.874130272238094), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (4.966656139420439, 1.161541590682067, 15.993038357924485, 7.175312328808405), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.21014591709747, 5.326768071513866, -18.450881488018126, 19.874779975491407), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (4.93668626355443, 1.1534174784892792, 15.99279010971033, 7.188033727750362), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.166644671858098, 5.277365365528824, -18.448718976425834, 19.851678563902254), 0.001 * arcsecond)

def test_mars_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 0.5912188976380217, 0.5 * meter)

    astrometric = e.observe(de405.mars)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 16.0296606272219, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.12731030858147, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.02988433983068, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.12820262180176, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.99950982315885, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.046277103674843, 0.001 * arcsecond)

def test_mars_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 1.4302506796029135, 0.5 * meter)

    astrometric = e.observe(de405.mars)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 23.545034875459514, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.882249043221036, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.544892038854186, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.882993630898111, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.536847630733252, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.935089760397493, 0.001 * arcsecond)

def test_mars_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 1.8496039270835372, 0.5 * meter)

    astrometric = e.observe(de405.mars)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 22.034936616343344, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.180707411034978, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.03468760932384, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.182134899635475, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.034417492807563, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.182689288940113, 0.001 * arcsecond)

def test_mars_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 1.766552316866877, 0.5 * meter)

    astrometric = e.observe(de405.mars)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 13.894324196598355, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.122808318928705, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.894132382683363, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.121796956140242, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.9057161859901, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.184654273116957, 0.001 * arcsecond)

def test_mars_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, (0.5912188976380217, 1.4302506796029135, 1.8496039270835372, 1.766552316866877), 0.5 * meter)

    astrometric = e.observe(de405.mars)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (16.0296606272219, 23.545034875459514, 22.034936616343344, 13.894324196598355), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-24.12731030858147, -4.882249043221036, -13.180707411034978, -12.122808318928705), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (16.02988433983068, 23.544892038854186, 22.03468760932384, 13.894132382683363), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-24.12820262180176, -4.882993630898111, -13.182134899635475, -12.121796956140242), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.99950982315885, 23.536847630733252, 22.034417492807563, 13.9057161859901), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-24.046277103674843, -4.935089760397493, -13.182689288940113, -12.184654273116957), 0.001 * arcsecond)

def test_jupiter_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 5.841600319231745, 0.5 * meter)

    astrometric = e.observe(de405.jupiter)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.104091505864655, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.6513409058207986, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.103936313614676, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.6524656208782568, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.07798204538282, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.8216129394812306, 0.001 * arcsecond)

def test_jupiter_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 5.913287883102949, 0.5 * meter)

    astrometric = e.observe(de405.jupiter)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 6.765154678701348, 0.001 * ra_arcsecond)
    compare(dec.degrees, 23.17039770012201, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 6.764854244708427, 0.001 * ra_arcsecond)
    compare(dec.degrees, 23.170736332068756, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 6.755383083025231, 0.001 * ra_arcsecond)
    compare(dec.degrees, 23.18268469367657, 0.001 * arcsecond)

def test_jupiter_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 4.621126565890217, 0.5 * meter)

    astrometric = e.observe(de405.jupiter)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 1.5913207023268698, 0.001 * ra_arcsecond)
    compare(dec.degrees, 8.595887646396902, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.5914167941833441, 0.001 * ra_arcsecond)
    compare(dec.degrees, 8.59631203599914, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.5911888424331277, 0.001 * ra_arcsecond)
    compare(dec.degrees, 8.594250857972389, 0.001 * arcsecond)

def test_jupiter_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 5.129958529243068, 0.5 * meter)

    astrometric = e.observe(de405.jupiter)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 4.822841055032964, 0.001 * ra_arcsecond)
    compare(dec.degrees, 21.649994488649472, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.822764769132737, 0.001 * ra_arcsecond)
    compare(dec.degrees, 21.64994169521301, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.835670404865467, 0.001 * ra_arcsecond)
    compare(dec.degrees, 21.670586389437943, 0.001 * arcsecond)

def test_jupiter_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, (5.841600319231745, 5.913287883102949, 4.621126565890217, 5.129958529243068), 0.5 * meter)

    astrometric = e.observe(de405.jupiter)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.104091505864655, 6.765154678701348, 1.5913207023268698, 4.822841055032964), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.6513409058207986, 23.17039770012201, 8.595887646396902, 21.649994488649472), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.103936313614676, 6.764854244708427, 1.5914167941833441, 4.822764769132737), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.6524656208782568, 23.170736332068756, 8.59631203599914, 21.64994169521301), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.07798204538282, 6.755383083025231, 1.5911888424331277, 4.835670404865467), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.8216129394812306, 23.18268469367657, 8.594250857972389, 21.670586389437943), 0.001 * arcsecond)

def test_saturn_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 9.382032444401027, 0.5 * meter)

    astrometric = e.observe(de405.saturn)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.462774885242021, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.045819985925936, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.4627075937035277, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.045735497802626, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.4352879582290172, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.911566107576899, 0.001 * arcsecond)

def test_saturn_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 9.420484451056101, 0.5 * meter)

    astrometric = e.observe(de405.saturn)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 19.814248756112033, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.93339019805076, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.814463444515557, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.932846451357467, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.805277718955743, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.958164640919687, 0.001 * arcsecond)

def test_saturn_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 8.652750126001484, 0.5 * meter)

    astrometric = e.observe(de405.saturn)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.584400980536592, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.616288735770388, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.584593321351076, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.616983167644806, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.5843611215084574, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.614774672730574, 0.001 * arcsecond)

def test_saturn_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 10.326368974662916, 0.5 * meter)

    astrometric = e.observe(de405.saturn)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 13.628484577191722, 0.001 * ra_arcsecond)
    compare(dec.degrees, -7.659435207931653, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.62827504244793, 0.001 * ra_arcsecond)
    compare(dec.degrees, -7.658028344724226, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.639628746850631, 0.001 * ra_arcsecond)
    compare(dec.degrees, -7.7232016421026195, 0.001 * arcsecond)

def test_saturn_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, (9.382032444401027, 9.420484451056101, 8.652750126001484, 10.326368974662916), 0.5 * meter)

    astrometric = e.observe(de405.saturn)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (2.462774885242021, 19.814248756112033, 2.584400980536592, 13.628484577191722), 0.001 * ra_arcsecond)
    compare(dec.degrees, (12.045819985925936, -20.93339019805076, 12.616288735770388, -7.659435207931653), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (2.4627075937035277, 19.814463444515557, 2.584593321351076, 13.62827504244793), 0.001 * ra_arcsecond)
    compare(dec.degrees, (12.045735497802626, -20.932846451357467, 12.616983167644806, -7.658028344724226), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (2.4352879582290172, 19.805277718955743, 2.5843611215084574, 13.639628746850631), 0.001 * ra_arcsecond)
    compare(dec.degrees, (11.911566107576899, -20.958164640919687, 12.614774672730574, -7.7232016421026195), 0.001 * arcsecond)

def test_uranus_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 18.75197906203834, 0.5 * meter)

    astrometric = e.observe(de405.uranus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.087167068351333, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.20723926118363262, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.087010426255668, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.20832526777272886, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.061052547705433, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.37749969290358604, 0.001 * arcsecond)

def test_uranus_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 18.62241700929518, 0.5 * meter)

    astrometric = e.observe(de405.uranus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.668551452013403, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.437331340689155, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.668859170516964, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.437016930580608, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.659361133085376, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.447681812488966, 0.001 * arcsecond)

def test_uranus_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 20.727159134679393, 0.5 * meter)

    astrometric = e.observe(de405.uranus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 21.16558686754142, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.018831731314233, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 21.165269485049027, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.020267168405788, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 21.164987614252272, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.020320613172, 0.001 * arcsecond)

def test_uranus_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 19.234768680195387, 0.5 * meter)

    astrometric = e.observe(de405.uranus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 0.48916431485643164, 0.001 * ra_arcsecond)
    compare(dec.degrees, 2.3565095329111823, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 0.48944632565389884, 0.001 * ra_arcsecond)
    compare(dec.degrees, 2.3583696385163115, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 0.5005500654503396, 0.001 * ra_arcsecond)
    compare(dec.degrees, 2.429779341040802, 0.001 * arcsecond)

def test_uranus_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, (18.75197906203834, 18.62241700929518, 20.727159134679393, 19.234768680195387), 0.5 * meter)

    astrometric = e.observe(de405.uranus)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.087167068351333, 18.668551452013403, 21.16558686754142, 0.48916431485643164), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.20723926118363262, -23.437331340689155, -17.018831731314233, 2.3565095329111823), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.087010426255668, 18.668859170516964, 21.165269485049027, 0.48944632565389884), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.20832526777272886, -23.437016930580608, -17.020267168405788, 2.3583696385163115), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.061052547705433, 18.659361133085376, 21.164987614252272, 0.5005500654503396), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.37749969290358604, -23.447681812488966, -17.020320613172, 2.429779341040802), 0.001 * arcsecond)

def test_neptune_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 29.83221264621946, 0.5 * meter)

    astrometric = e.observe(de405.neptune)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 15.637210587139663, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.67999613660563, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.63739098768298, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.680453730264624, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.608486730597074, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.58379328551931, 0.001 * arcsecond)

def test_neptune_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 29.490001740438892, 0.5 * meter)

    astrometric = e.observe(de405.neptune)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 19.03623522579387, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.792864018500975, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.036513633320563, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.792510662370386, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.02716408230529, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.808047913986798, 0.001 * arcsecond)

def test_neptune_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 31.02449192035449, 0.5 * meter)

    astrometric = e.observe(de405.neptune)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 20.362841834121518, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.212425239376326, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 20.36247543901059, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.213645950878377, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 20.36218815756048, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.213233798897665, 0.001 * arcsecond)

def test_neptune_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 28.98411802971635, 0.5 * meter)

    astrometric = e.observe(de405.neptune)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 22.252468120719442, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.504657215501586, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.252825961036415, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.502649482645886, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.2643158309744, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.437330191299896, 0.001 * arcsecond)

def test_neptune_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, (29.83221264621946, 29.490001740438892, 31.02449192035449, 28.98411802971635), 0.5 * meter)

    astrometric = e.observe(de405.neptune)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (15.637210587139663, 19.03623522579387, 20.362841834121518, 22.252468120719442), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-17.67999613660563, -21.792864018500975, -19.212425239376326, -11.504657215501586), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (15.63739098768298, 19.036513633320563, 20.36247543901059, 22.252825961036415), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-17.680453730264624, -21.792510662370386, -19.213645950878377, -11.502649482645886), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.608486730597074, 19.02716408230529, 20.36218815756048, 22.2643158309744), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-17.58379328551931, -21.808047913986798, -19.213233798897665, -11.437330191299896), 0.001 * arcsecond)

def test_pluto_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 32.312971776632494, 0.5 * meter)

    astrometric = e.observe(de405.pluto)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.015311208821212, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.620557180992584, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.01514128380381, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.62299016066861, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 11.989232654068259, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.792242650891883, 0.001 * arcsecond)

def test_pluto_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 28.707485955458118, 0.5 * meter)

    astrometric = e.observe(de405.pluto)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 15.216302246424346, 0.001 * ra_arcsecond)
    compare(dec.degrees, -1.3346560528819575, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.216661036271791, 0.001 * ra_arcsecond)
    compare(dec.degrees, -1.3358630622052712, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.208581663980876, 0.001 * ra_arcsecond)
    compare(dec.degrees, -1.3022394883151638, 0.001 * arcsecond)

def test_pluto_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 31.064412196006607, 0.5 * meter)

    astrometric = e.observe(de405.pluto)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 16.761873062250743, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.39643313463007, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.76152667540677, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.396301545071509, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.761277438459963, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.394288734411235, 0.001 * arcsecond)

def test_pluto_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 31.69909782133192, 0.5 * meter)

    astrometric = e.observe(de405.pluto)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.488351288595236, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.552190994888846, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.488573622605898, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.551729414764313, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.501338273669152, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.54122790974374, 0.001 * arcsecond)

def test_pluto_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, (32.312971776632494, 28.707485955458118, 31.064412196006607, 31.69909782133192), 0.5 * meter)

    astrometric = e.observe(de405.pluto)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.015311208821212, 15.216302246424346, 16.761873062250743, 18.488351288595236), 0.001 * ra_arcsecond)
    compare(dec.degrees, (16.620557180992584, -1.3346560528819575, -11.39643313463007, -19.552190994888846), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.01514128380381, 15.216661036271791, 16.76152667540677, 18.488573622605898), 0.001 * ra_arcsecond)
    compare(dec.degrees, (16.62299016066861, -1.3358630622052712, -11.396301545071509, -19.551729414764313), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (11.989232654068259, 15.208581663980876, 16.761277438459963, 18.501338273669152), 0.001 * ra_arcsecond)
    compare(dec.degrees, (16.792242650891883, -1.3022394883151638, -11.394288734411235, -19.54122790974374), 0.001 * arcsecond)

def test_sun_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 1.0160878650466754, 0.5 * meter)

    astrometric = e.observe(de405.sun)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 8.03008088792976, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.496475643233936, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 8.029690303049978, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.49760546326072, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 8.000108116572395, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.584930935996038, 0.001 * arcsecond)

def test_sun_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 1.0118605934887042, 0.5 * meter)

    astrometric = e.observe(de405.sun)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 3.776110727862678, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.907832379364578, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 3.7757213854872145, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.906601181542, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 3.7666292045824337, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.87917377230974, 0.001 * arcsecond)

def test_sun_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 0.9833276788862821, 0.5 * meter)

    astrometric = e.observe(de405.sun)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 18.752544254682526, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.033309607967187, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.752126228091367, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.03376015263556, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.75183797477899, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.032488638722818, 0.001 * arcsecond)

def test_sun_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 1.0107820040799866, 0.5 * meter)

    astrometric = e.observe(de405.sun)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 10.268162490439074, 0.001 * ra_arcsecond)
    compare(dec.degrees, 10.751933902906117, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 10.267805651450436, 0.001 * ra_arcsecond)
    compare(dec.degrees, 10.753946960547601, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 10.279264504672037, 0.001 * ra_arcsecond)
    compare(dec.degrees, 10.68850786534133, 0.001 * arcsecond)

def test_sun_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, (1.0160878650466754, 1.0118605934887042, 0.9833276788862821, 1.0107820040799866), 0.5 * meter)

    astrometric = e.observe(de405.sun)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (8.03008088792976, 3.776110727862678, 18.752544254682526, 10.268162490439074), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.496475643233936, 19.907832379364578, -23.033309607967187, 10.751933902906117), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (8.029690303049978, 3.7757213854872145, 18.752126228091367, 10.267805651450436), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.49760546326072, 19.906601181542, -23.03376015263556, 10.753946960547601), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (8.000108116572395, 3.7666292045824337, 18.75183797477899, 10.279264504672037), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.584930935996038, 19.87917377230974, -23.032488638722818, 10.68850786534133), 0.001 * arcsecond)

def test_moon_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)

    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.0026034424248855136, 0.5 * meter)

    astrometric = e.observe(de405.moon)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 12.472463241145173, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.546618838171283, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.47234028706646, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.545964408924448, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.446262111681095, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.378227942513377, 0.001 * arcsecond)

def test_moon_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.002481509229659986, 0.5 * meter)

    astrometric = e.observe(de405.moon)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 23.67644381740948, 0.001 * ra_arcsecond)
    compare(dec.degrees, 1.8587554901327863, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.676289920709504, 0.001 * ra_arcsecond)
    compare(dec.degrees, 1.8574138759902254, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.668278096873856, 0.001 * ra_arcsecond)
    compare(dec.degrees, 1.805189185726724, 0.001 * arcsecond)

def test_moon_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)

    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.0026902029885132865, 0.5 * meter)

    astrometric = e.observe(de405.moon)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 14.830020573942209, 0.001 * ra_arcsecond)
    compare(dec.degrees, -10.900635500944452, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 14.829807890359646, 0.001 * ra_arcsecond)
    compare(dec.degrees, -10.900127758842366, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 14.82957327176072, 0.001 * ra_arcsecond)
    compare(dec.degrees, -10.897905576905867, 0.001 * arcsecond)

def test_moon_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)

    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.0024739078649309238, 0.5 * meter)

    astrometric = e.observe(de405.moon)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 16.39102815233177, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.93676001523414, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.39106196861365, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.936774891979844, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.403831131432188, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.96508913558473, 0.001 * arcsecond)

def test_moon_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)

    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, (0.0026034424248855136, 0.002481509229659986, 0.0026902029885132865, 0.0024739078649309238), 0.5 * meter)

    astrometric = e.observe(de405.moon)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (12.472463241145173, 23.67644381740948, 14.830020573942209, 16.39102815233177), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-4.546618838171283, 1.8587554901327863, -10.900635500944452, -20.93676001523414), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.47234028706646, 23.676289920709504, 14.829807890359646, 16.39106196861365), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-4.545964408924448, 1.8574138759902254, -10.900127758842366, -20.936774891979844), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.446262111681095, 23.668278096873856, 14.82957327176072, 16.403831131432188), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-4.378227942513377, 1.805189185726724, -10.897905576905867, -20.96508913558473), 0.001 * arcsecond)

def test_polaris_geocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    e = de405.earth(jd)
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.528369749952934, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26420848455291, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.522801492978089, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.25882879505869, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.0385816433557027, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.11999387030946, 0.001 * arcsecond)

def test_polaris_geocentric_date1():
    jd = JulianDate(tt=2448031.5)
    e = de405.earth(jd)
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.529691027594406, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26413894692217, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.5033568528110783, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26201007627152, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.332921180528844, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.22082922133737, 0.001 * arcsecond)

def test_polaris_geocentric_date2():
    jd = JulianDate(tt=2451545.0)
    e = de405.earth(jd)
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.5302921882000122, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26411027119273, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.5446332154627274, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26917874902797, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.5459982729094506, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26697328449004, 0.001 * arcsecond)

def test_polaris_geocentric_date3():
    jd = JulianDate(tt=2456164.5)
    e = de405.earth(jd)
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, 2.5311170656101485, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26406906493732, 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.5416095337355356, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.25923373182653, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.8064741334456444, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.3136939266471, 0.001 * arcsecond)

def test_polaris_geocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    e = de405.earth(jd)
    star = starlib.Star(ra_hours=2.530301028, dec_degrees=89.264109444,
                        ra_mas_per_year=44.22, dec_mas_per_year=-11.75,
                        parallax_mas=7.56, radial_km_per_s=-17.4)

    astrometric = e.observe(star)
    ra, dec, distance = astrometric.radec()
    compare(ra.hours, (2.528369749952934, 2.529691027594406, 2.5302921882000122, 2.5311170656101485), 0.001 * ra_arcsecond)
    compare(dec.degrees, (89.26420848455291, 89.26413894692217, 89.26411027119273, 89.26406906493732), 0.001 * arcsecond)

    apparent = astrometric.apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (2.522801492978089, 2.5033568528110783, 2.5446332154627274, 2.5416095337355356), 0.001 * ra_arcsecond)
    compare(dec.degrees, (89.25882879505869, 89.26201007627152, 89.26917874902797, 89.25923373182653), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (2.0385816433557027, 2.332921180528844, 2.5459982729094506, 2.8064741334456444), 0.001 * ra_arcsecond)
    compare(dec.degrees, (89.11999387030946, 89.22082922133737, 89.26697328449004, 89.3136939266471), 0.001 * arcsecond)

def test_mercury_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mercury).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.9049140222444105, 0.001 * ra_arcsecond)
    compare(dec.degrees, 22.33276016366846, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.874898511438326, 0.001 * ra_arcsecond)
    compare(dec.degrees, 22.41529463722477, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 46.32122675660419, 0.001 * arcsecond)
    compare(az.degrees, 262.1859052156761, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 46.33688339908464, 0.001 * arcsecond)
    compare(az.degrees, 262.1859052156761, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 46.33704240111, 0.001 * arcsecond)
    compare(az.degrees, 262.1859052156761, 0.001 * arcsecond)

def test_mercury_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mercury).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.4699595920648565, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.24594905426479, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.461508188066483, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.205182598299665, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -17.340667089884377, 0.001 * arcsecond)
    compare(az.degrees, 300.9176579181716, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -17.340667089884377, 0.001 * arcsecond)
    compare(az.degrees, 300.9176579181716, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -17.340667089884377, 0.001 * arcsecond)
    compare(az.degrees, 300.9176579181716, 0.001 * arcsecond)

def test_mercury_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mercury).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.138603904058247, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.421550562485436, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.138318996641566, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.420003066967507, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -0.12765060376707993, 0.001 * arcsecond)
    compare(az.degrees, 121.97764361867154, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 0.36890915770102595, 0.001 * arcsecond)
    compare(az.degrees, 121.97764361867154, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 0.3731892291678207, 0.001 * arcsecond)
    compare(az.degrees, 121.97764361867154, 0.001 * arcsecond)

def test_mercury_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mercury).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 9.295468142561822, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.685908124650233, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 9.307459135231527, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.630243128506475, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -9.116616855755936, 0.001 * arcsecond)
    compare(az.degrees, 300.14202643731034, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -9.116616855755936, 0.001 * arcsecond)
    compare(az.degrees, 300.14202643731034, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -9.116616855755936, 0.001 * arcsecond)
    compare(az.degrees, 300.14202643731034, 0.001 * arcsecond)

def test_mercury_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mercury).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (7.9049140222444105, 2.4699595920648565, 18.138603904058247, 9.295468142561822), 0.001 * ra_arcsecond)
    compare(dec.degrees, (22.33276016366846, 11.24594905426479, -24.421550562485436, 16.685908124650233), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (7.874898511438326, 2.461508188066483, 18.138318996641566, 9.307459135231527), 0.001 * ra_arcsecond)
    compare(dec.degrees, (22.41529463722477, 11.205182598299665, -24.420003066967507, 16.630243128506475), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (46.32122675660419, -17.340667089884377, -0.12765060376707993, -9.116616855755936), 0.001 * arcsecond)
    compare(az.degrees, (262.1859052156761, 300.9176579181716, 121.97764361867154, 300.14202643731034), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (46.33688339908464, -17.340667089884377, 0.36890915770102595, -9.116616855755936), 0.001 * arcsecond)
    compare(az.degrees, (262.1859052156761, 300.9176579181716, 121.97764361867154, 300.14202643731034), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (46.33704240111, -17.340667089884377, 0.3731892291678207, -9.116616855755936), 0.001 * arcsecond)
    compare(az.degrees, (262.1859052156761, 300.9176579181716, 121.97764361867154, 300.14202643731034), 0.001 * arcsecond)

def test_venus_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.venus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.9665155792599744, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.208668727034965, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.936546062416393, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.165161469755116, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 11.152374062991527, 0.001 * arcsecond)
    compare(az.degrees, 287.0030740239525, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 11.231992752470703, 0.001 * arcsecond)
    compare(az.degrees, 287.0030740239525, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 11.23279626216302, 0.001 * arcsecond)
    compare(az.degrees, 287.0030740239525, 0.001 * arcsecond)

def test_venus_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.venus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.1614662937271139, 0.001 * ra_arcsecond)
    compare(dec.degrees, 5.3252225859555455, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.1533422187037872, 0.001 * ra_arcsecond)
    compare(dec.degrees, 5.275819541572406, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -34.134914076462266, 0.001 * arcsecond)
    compare(az.degrees, 313.64872862118426, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -34.134914076462266, 0.001 * arcsecond)
    compare(az.degrees, 313.64872862118426, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -34.134914076462266, 0.001 * arcsecond)
    compare(az.degrees, 313.64872862118426, 0.001 * arcsecond)

def test_venus_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.venus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.99311221167692, 0.001 * ra_arcsecond)
    compare(dec.degrees, -18.452566802886178, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.992863961375889, 0.001 * ra_arcsecond)
    compare(dec.degrees, -18.450404301558013, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 23.228910604670872, 0.001 * arcsecond)
    compare(az.degrees, 142.11613981416264, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 23.266773672986062, 0.001 * arcsecond)
    compare(az.degrees, 142.11613981416264, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 23.267157712313733, 0.001 * arcsecond)
    compare(az.degrees, 142.11613981416264, 0.001 * arcsecond)

def test_venus_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.venus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 7.175218975921812, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.87224931182422, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 7.187940160922054, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.84914957337174, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -24.35999541091543, 0.001 * arcsecond)
    compare(az.degrees, 327.640588969984, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -24.35999541091543, 0.001 * arcsecond)
    compare(az.degrees, 327.640588969984, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -24.35999541091543, 0.001 * arcsecond)
    compare(az.degrees, 327.640588969984, 0.001 * arcsecond)

def test_venus_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.venus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (4.9665155792599744, 1.1614662937271139, 15.99311221167692, 7.175218975921812), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.208668727034965, 5.3252225859555455, -18.452566802886178, 19.87224931182422), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (4.936546062416393, 1.1533422187037872, 15.992863961375889, 7.187940160922054), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.165161469755116, 5.275819541572406, -18.450404301558013, 19.84914957337174), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (11.152374062991527, -34.134914076462266, 23.228910604670872, -24.35999541091543), 0.001 * arcsecond)
    compare(az.degrees, (287.0030740239525, 313.64872862118426, 142.11613981416264, 327.640588969984), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (11.231992752470703, -34.134914076462266, 23.266773672986062, -24.35999541091543), 0.001 * arcsecond)
    compare(az.degrees, (287.0030740239525, 313.64872862118426, 142.11613981416264, 327.640588969984), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (11.23279626216302, -34.134914076462266, 23.267157712313733, -24.35999541091543), 0.001 * arcsecond)
    compare(az.degrees, (287.0030740239525, 313.64872862118426, 142.11613981416264, 327.640588969984), 0.001 * arcsecond)

def test_mars_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mars).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.030112454663165, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.13088318769704, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.999737237126766, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.04896650222992, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -3.540294697029495, 0.001 * arcsecond)
    compare(az.degrees, 118.34877634707443, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -3.540294697029495, 0.001 * arcsecond)
    compare(az.degrees, 118.34877634707443, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -3.540294697029495, 0.001 * arcsecond)
    compare(az.degrees, 118.34877634707443, 0.001 * arcsecond)

def test_mars_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mars).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.54486790147113, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.883946644223004, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.53682348628842, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.93604274443558, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -54.1089628741949, 0.001 * arcsecond)
    compare(az.degrees, 338.0117138951488, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -54.1089628741949, 0.001 * arcsecond)
    compare(az.degrees, 338.0117138951488, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -54.1089628741949, 0.001 * arcsecond)
    compare(az.degrees, 338.0117138951488, 0.001 * arcsecond)

def test_mars_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mars).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.034740913364253, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.182784253332379, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.03447079524992, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.183338672731743, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -36.90573266459917, 0.001 * arcsecond)
    compare(az.degrees, 76.12368450672824, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -36.90573266459917, 0.001 * arcsecond)
    compare(az.degrees, 76.12368450672824, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -36.90573266459917, 0.001 * arcsecond)
    compare(az.degrees, 76.12368450672824, 0.001 * arcsecond)

def test_mars_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mars).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.8940809044733, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.122804110106658, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.905664739133574, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.185661905051246, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 22.094794272017666, 0.001 * arcsecond)
    compare(az.degrees, 231.6381663847761, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 22.134776069489533, 0.001 * arcsecond)
    compare(az.degrees, 231.6381663847761, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 22.135181528743814, 0.001 * arcsecond)
    compare(az.degrees, 231.6381663847761, 0.001 * arcsecond)

def test_mars_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.mars).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (16.030112454663165, 23.54486790147113, 22.034740913364253, 13.8940809044733), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-24.13088318769704, -4.883946644223004, -13.182784253332379, -12.122804110106658), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.999737237126766, 23.53682348628842, 22.03447079524992, 13.905664739133574), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-24.04896650222992, -4.93604274443558, -13.183338672731743, -12.185661905051246), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (-3.540294697029495, -54.1089628741949, -36.90573266459917, 22.094794272017666), 0.001 * arcsecond)
    compare(az.degrees, (118.34877634707443, 338.0117138951488, 76.12368450672824, 231.6381663847761), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (-3.540294697029495, -54.1089628741949, -36.90573266459917, 22.134776069489533), 0.001 * arcsecond)
    compare(az.degrees, (118.34877634707443, 338.0117138951488, 76.12368450672824, 231.6381663847761), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (-3.540294697029495, -54.1089628741949, -36.90573266459917, 22.135181528743814), 0.001 * arcsecond)
    compare(az.degrees, (118.34877634707443, 338.0117138951488, 76.12368450672824, 231.6381663847761), 0.001 * arcsecond)

def test_jupiter_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.jupiter).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.103946503374882, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.6522085918269473, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.077992233588102, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.8213558931137472, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 49.406516031446415, 0.001 * arcsecond)
    compare(az.degrees, 156.07088561561812, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 49.42056980196561, 0.001 * arcsecond)
    compare(az.degrees, 156.07088561561812, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 49.420712533159296, 0.001 * arcsecond)
    compare(az.degrees, 156.07088561561812, 0.001 * arcsecond)

def test_jupiter_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.jupiter).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 6.76483682133995, 0.001 * ra_arcsecond)
    compare(dec.degrees, 23.170587900559504, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 6.755365668515656, 0.001 * ra_arcsecond)
    compare(dec.degrees, 23.182536029964222, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 38.005051266909966, 0.001 * arcsecond)
    compare(az.degrees, 270.63795554820535, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 38.02600464378365, 0.001 * arcsecond)
    compare(az.degrees, 270.63795554820535, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 38.026217393249304, 0.001 * arcsecond)
    compare(az.degrees, 270.63795554820535, 0.001 * arcsecond)

def test_jupiter_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.jupiter).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 1.5914118935512864, 0.001 * ra_arcsecond)
    compare(dec.degrees, 8.595923929888198, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 1.5911839414385693, 0.001 * ra_arcsecond)
    compare(dec.degrees, 8.593862752942394, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -42.482560972481394, 0.001 * arcsecond)
    compare(az.degrees, 359.3596746827537, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -42.482560972481394, 0.001 * arcsecond)
    compare(az.degrees, 359.3596746827537, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -42.482560972481394, 0.001 * arcsecond)
    compare(az.degrees, 359.3596746827537, 0.001 * arcsecond)

def test_jupiter_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.jupiter).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 4.82276173655752, 0.001 * ra_arcsecond)
    compare(dec.degrees, 21.649526689253502, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 4.835667333191383, 0.001 * ra_arcsecond)
    compare(dec.degrees, 21.670171438742262, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -29.289013841967986, 0.001 * arcsecond)
    compare(az.degrees, 4.327425566855531, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -29.289013841967986, 0.001 * arcsecond)
    compare(az.degrees, 4.327425566855531, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -29.289013841967986, 0.001 * arcsecond)
    compare(az.degrees, 4.327425566855531, 0.001 * arcsecond)

def test_jupiter_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.jupiter).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.103946503374882, 6.76483682133995, 1.5914118935512864, 4.82276173655752), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.6522085918269473, 23.170587900559504, 8.595923929888198, 21.649526689253502), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.077992233588102, 6.755365668515656, 1.5911839414385693, 4.835667333191383), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.8213558931137472, 23.182536029964222, 8.593862752942394, 21.670171438742262), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (49.406516031446415, 38.005051266909966, -42.482560972481394, -29.289013841967986), 0.001 * arcsecond)
    compare(az.degrees, (156.07088561561812, 270.63795554820535, 359.3596746827537, 4.327425566855531), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (49.42056980196561, 38.02600464378365, -42.482560972481394, -29.289013841967986), 0.001 * arcsecond)
    compare(az.degrees, (156.07088561561812, 270.63795554820535, 359.3596746827537, 4.327425566855531), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (49.420712533159296, 38.026217393249304, -42.482560972481394, -29.289013841967986), 0.001 * arcsecond)
    compare(az.degrees, (156.07088561561812, 270.63795554820535, 359.3596746827537, 4.327425566855531), 0.001 * arcsecond)

def test_saturn_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.saturn).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.46269388589056, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.045561201575385, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.435274279115236, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.911391441362449, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -20.66268694032327, 0.001 * arcsecond)
    compare(az.degrees, 306.01978569992684, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -20.66268694032327, 0.001 * arcsecond)
    compare(az.degrees, 306.01978569992684, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -20.66268694032327, 0.001 * arcsecond)
    compare(az.degrees, 306.01978569992684, 0.001 * arcsecond)

def test_saturn_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.saturn).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.814469727768646, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.932928080758657, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.805283998285297, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.958246345579152, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -48.93337647838982, 0.001 * arcsecond)
    compare(az.degrees, 76.8837444919445, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -48.93337647838982, 0.001 * arcsecond)
    compare(az.degrees, 76.8837444919445, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -48.93337647838982, 0.001 * arcsecond)
    compare(az.degrees, 76.8837444919445, 0.001 * arcsecond)

def test_saturn_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.saturn).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 2.5845847757319116, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.616768688416157, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 2.5843525758885217, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.614560194137898, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -36.501918751911674, 0.001 * arcsecond)
    compare(az.degrees, 341.22347230453323, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -36.501918751911674, 0.001 * arcsecond)
    compare(az.degrees, 341.22347230453323, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -36.501918751911674, 0.001 * arcsecond)
    compare(az.degrees, 341.22347230453323, 0.001 * arcsecond)

def test_saturn_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.saturn).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 13.628268137367913, 0.001 * ra_arcsecond)
    compare(dec.degrees, -7.658197329820582, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 13.639621846921335, 0.001 * ra_arcsecond)
    compare(dec.degrees, -7.723370683249702, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 22.966758516111867, 0.001 * arcsecond)
    compare(az.degrees, 238.00627672875675, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 23.005094362956058, 0.001 * arcsecond)
    compare(az.degrees, 238.00627672875675, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 23.005483182929083, 0.001 * arcsecond)
    compare(az.degrees, 238.00627672875675, 0.001 * arcsecond)

def test_saturn_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.saturn).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (2.46269388589056, 19.814469727768646, 2.5845847757319116, 13.628268137367913), 0.001 * ra_arcsecond)
    compare(dec.degrees, (12.045561201575385, -20.932928080758657, 12.616768688416157, -7.658197329820582), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (2.435274279115236, 19.805283998285297, 2.5843525758885217, 13.639621846921335), 0.001 * ra_arcsecond)
    compare(dec.degrees, (11.911391441362449, -20.958246345579152, 12.614560194137898, -7.723370683249702), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (-20.66268694032327, -48.93337647838982, -36.501918751911674, 22.966758516111867), 0.001 * arcsecond)
    compare(az.degrees, (306.01978569992684, 76.8837444919445, 341.22347230453323, 238.00627672875675), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (-20.66268694032327, -48.93337647838982, -36.501918751911674, 23.005094362956058), 0.001 * arcsecond)
    compare(az.degrees, (306.01978569992684, 76.8837444919445, 341.22347230453323, 238.00627672875675), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (-20.66268694032327, -48.93337647838982, -36.501918751911674, 23.005483182929083), 0.001 * arcsecond)
    compare(az.degrees, (306.01978569992684, 76.8837444919445, 341.22347230453323, 238.00627672875675), 0.001 * arcsecond)

def test_uranus_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.uranus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.087016642067397, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.2082444210471117, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.061058763070791, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.377418836834601, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 49.06396822144691, 0.001 * arcsecond)
    compare(az.degrees, 156.6525604020511, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 49.07819253506017, 0.001 * arcsecond)
    compare(az.degrees, 156.6525604020511, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 49.078336997561024, 0.001 * arcsecond)
    compare(az.degrees, 156.6525604020511, 0.001 * arcsecond)

def test_uranus_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.uranus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.668863148648313, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.43704804377174, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.65936510933368, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.447712978993902, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -37.0259637798912, 0.001 * arcsecond)
    compare(az.degrees, 91.80748703145905, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -37.0259637798912, 0.001 * arcsecond)
    compare(az.degrees, 91.80748703145905, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -37.0259637798912, 0.001 * arcsecond)
    compare(az.degrees, 91.80748703145905, 0.001 * arcsecond)

def test_uranus_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.uranus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 21.16527335872666, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.02030811911839, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 21.164991487815, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.02036156614209, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -29.175475562665554, 0.001 * arcsecond)
    compare(az.degrees, 88.85671230431439, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -29.175475562665554, 0.001 * arcsecond)
    compare(az.degrees, 88.85671230431439, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -29.175475562665554, 0.001 * arcsecond)
    compare(az.degrees, 88.85671230431439, 0.001 * arcsecond)

def test_uranus_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.uranus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 0.48945083888242813, 0.001 * ra_arcsecond)
    compare(dec.degrees, 2.3582861967255484, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 0.5005545778925, 0.001 * ra_arcsecond)
    compare(dec.degrees, 2.429695886841979, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -14.5260443119261, 0.001 * arcsecond)
    compare(az.degrees, 74.60219420538263, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -14.5260443119261, 0.001 * arcsecond)
    compare(az.degrees, 74.60219420538263, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -14.5260443119261, 0.001 * arcsecond)
    compare(az.degrees, 74.60219420538263, 0.001 * arcsecond)

def test_uranus_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.uranus).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.087016642067397, 18.668863148648313, 21.16527335872666, 0.48945083888242813), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.2082444210471117, -23.43704804377174, -17.02030811911839, 2.3582861967255484), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.061058763070791, 18.65936510933368, 21.164991487815, 0.5005545778925), 0.001 * ra_arcsecond)
    compare(dec.degrees, (0.377418836834601, -23.447712978993902, -17.02036156614209, 2.429695886841979), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (49.06396822144691, -37.0259637798912, -29.175475562665554, -14.5260443119261), 0.001 * arcsecond)
    compare(az.degrees, (156.6525604020511, 91.80748703145905, 88.85671230431439, 74.60219420538263), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (49.07819253506017, -37.0259637798912, -29.175475562665554, -14.5260443119261), 0.001 * arcsecond)
    compare(az.degrees, (156.6525604020511, 91.80748703145905, 88.85671230431439, 74.60219420538263), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (49.078336997561024, -37.0259637798912, -29.175475562665554, -14.5260443119261), 0.001 * arcsecond)
    compare(az.degrees, (156.6525604020511, 91.80748703145905, 88.85671230431439, 74.60219420538263), 0.001 * arcsecond)

def test_neptune_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.neptune).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.637396931781986, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.6804899511715, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.60849266504413, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.583829722494016, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 4.869377826364499, 0.001 * arcsecond)
    compare(az.degrees, 117.29043762875322, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 5.031511017144567, 0.001 * arcsecond)
    compare(az.degrees, 117.29043762875322, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 5.033116634142289, 0.001 * arcsecond)
    compare(az.degrees, 117.29043762875322, 0.001 * arcsecond)

def test_neptune_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.neptune).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 19.036514568239326, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.792523874854826, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 19.027165016434417, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.808061138689617, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -40.43318694811052, 0.001 * arcsecond)
    compare(az.degrees, 86.51833613444356, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -40.43318694811052, 0.001 * arcsecond)
    compare(az.degrees, 86.51833613444356, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -40.43318694811052, 0.001 * arcsecond)
    compare(az.degrees, 86.51833613444356, 0.001 * arcsecond)

def test_neptune_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.neptune).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 20.362478654099593, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.21366591391132, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 20.36219137258442, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.21325376377245, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -21.10215467278759, 0.001 * arcsecond)
    compare(az.degrees, 98.1496208151544, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -21.10215467278759, 0.001 * arcsecond)
    compare(az.degrees, 98.1496208151544, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -21.10215467278759, 0.001 * arcsecond)
    compare(az.degrees, 98.1496208151544, 0.001 * arcsecond)

def test_neptune_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.neptune).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 22.252831344843074, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.502690543226889, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 22.264321215062385, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.4373712085964, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 2.4167829049998346, 0.001 * arcsecond)
    compare(az.degrees, 106.80925972576064, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 2.6713913487619294, 0.001 * arcsecond)
    compare(az.degrees, 106.80925972576064, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 2.6738334093304843, 0.001 * arcsecond)
    compare(az.degrees, 106.80925972576064, 0.001 * arcsecond)

def test_neptune_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.neptune).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (15.637396931781986, 19.036514568239326, 20.362478654099593, 22.252831344843074), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-17.6804899511715, -21.792523874854826, -19.21366591391132, -11.502690543226889), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (15.60849266504413, 19.027165016434417, 20.36219137258442, 22.264321215062385), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-17.583829722494016, -21.808061138689617, -19.21325376377245, -11.4373712085964), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (4.869377826364499, -40.43318694811052, -21.10215467278759, 2.4167829049998346), 0.001 * arcsecond)
    compare(az.degrees, (117.29043762875322, 86.51833613444356, 98.1496208151544, 106.80925972576064), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (5.031511017144567, -40.43318694811052, -21.10215467278759, 2.6713913487619294), 0.001 * arcsecond)
    compare(az.degrees, (117.29043762875322, 86.51833613444356, 98.1496208151544, 106.80925972576064), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (5.033116634142289, -40.43318694811052, -21.10215467278759, 2.6738334093304843), 0.001 * arcsecond)
    compare(az.degrees, (117.29043762875322, 86.51833613444356, 98.1496208151544, 106.80925972576064), 0.001 * arcsecond)

def test_pluto_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.pluto).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.015146948702718, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.622956629676775, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 11.989238323883423, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.792209116103166, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 64.72856074651932, 0.001 * arcsecond)
    compare(az.degrees, 147.2138070056032, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 64.73630449169255, 0.001 * arcsecond)
    compare(az.degrees, 147.2138070056032, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 64.7363831493004, 0.001 * arcsecond)
    compare(az.degrees, 147.2138070056032, 0.001 * arcsecond)

def test_pluto_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.pluto).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 15.216666873470118, 0.001 * ra_arcsecond)
    compare(dec.degrees, -1.335915234746897, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 15.208587498665663, 0.001 * ra_arcsecond)
    compare(dec.degrees, -1.30229172206482, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 16.2337344521234, 0.001 * arcsecond)
    compare(az.degrees, 105.39943656311958, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 16.288892801912894, 0.001 * arcsecond)
    compare(az.degrees, 105.39943656311958, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 16.28945132964904, 0.001 * arcsecond)
    compare(az.degrees, 105.39943656311958, 0.001 * arcsecond)

def test_pluto_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.pluto).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.761532920101487, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.396347593297175, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.76128368305737, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.394334784193745, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 22.700996363632996, 0.001 * arcsecond)
    compare(az.degrees, 127.81134408260581, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 22.739821647292274, 0.001 * arcsecond)
    compare(az.degrees, 127.81134408260581, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 22.74021541578692, 0.001 * arcsecond)
    compare(az.degrees, 127.81134408260581, 0.001 * arcsecond)

def test_pluto_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.pluto).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.488579709427018, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.551785355075804, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.501344365322606, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.541283736216645, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 28.339822068789154, 0.001 * arcsecond)
    compare(az.degrees, 157.51785266272373, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 28.37007124206125, 0.001 * arcsecond)
    compare(az.degrees, 157.51785266272373, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 28.370378222043676, 0.001 * arcsecond)
    compare(az.degrees, 157.51785266272373, 0.001 * arcsecond)

def test_pluto_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.pluto).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.015146948702718, 15.216666873470118, 16.761532920101487, 18.488579709427018), 0.001 * ra_arcsecond)
    compare(dec.degrees, (16.622956629676775, -1.335915234746897, -11.396347593297175, -19.551785355075804), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (11.989238323883423, 15.208587498665663, 16.76128368305737, 18.501344365322606), 0.001 * ra_arcsecond)
    compare(dec.degrees, (16.792209116103166, -1.30229172206482, -11.394334784193745, -19.541283736216645), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (64.72856074651932, 16.2337344521234, 22.700996363632996, 28.339822068789154), 0.001 * arcsecond)
    compare(az.degrees, (147.2138070056032, 105.39943656311958, 127.81134408260581, 157.51785266272373), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (64.73630449169255, 16.288892801912894, 22.739821647292274, 28.37007124206125), 0.001 * arcsecond)
    compare(az.degrees, (147.2138070056032, 105.39943656311958, 127.81134408260581, 157.51785266272373), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (64.7363831493004, 16.28945132964904, 22.74021541578692, 28.370378222043676), 0.001 * arcsecond)
    compare(az.degrees, (147.2138070056032, 105.39943656311958, 127.81134408260581, 157.51785266272373), 0.001 * arcsecond)

def test_sun_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.sun).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 8.02959789881544, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.496678572125123, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 8.000015838288707, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.584000539289498, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 46.72403357148921, 0.001 * arcsecond)
    compare(az.degrees, 258.5550717845947, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 46.73947196634785, 0.001 * arcsecond)
    compare(az.degrees, 258.5550717845947, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 46.73962875307822, 0.001 * arcsecond)
    compare(az.degrees, 258.5550717845947, 0.001 * arcsecond)

def test_sun_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.sun).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 3.775590638161117, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.90505409109931, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 3.76649857059908, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.87762515818775, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 2.2209469369832675, 0.001 * arcsecond)
    compare(az.degrees, 293.95636637272145, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 2.486840978779398, 0.001 * arcsecond)
    compare(az.degrees, 293.95636637272145, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 2.489379891081043, 0.001 * arcsecond)
    compare(az.degrees, 293.95636637272145, 0.001 * arcsecond)

def test_sun_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.sun).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 18.752264357691004, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.03532101826747, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 18.751976099155204, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.034049570458144, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -5.486505415022833, 0.001 * arcsecond)
    compare(az.degrees, 115.32008451470388, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -5.486505415022833, 0.001 * arcsecond)
    compare(az.degrees, 115.32008451470388, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -5.486505415022833, 0.001 * arcsecond)
    compare(az.degrees, 115.32008451470388, 0.001 * arcsecond)

def test_sun_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.sun).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 10.267679924967121, 0.001 * ra_arcsecond)
    compare(dec.degrees, 10.75239953710826, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 10.279138748598198, 0.001 * ra_arcsecond)
    compare(dec.degrees, 10.68696144441038, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -2.738407691502772, 0.001 * arcsecond)
    compare(az.degrees, 286.09632001391725, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -2.738407691502772, 0.001 * arcsecond)
    compare(az.degrees, 286.09632001391725, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -2.738407691502772, 0.001 * arcsecond)
    compare(az.degrees, 286.09632001391725, 0.001 * arcsecond)

def test_sun_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.sun).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (8.02959789881544, 3.775590638161117, 18.752264357691004, 10.267679924967121), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.496678572125123, 19.90505409109931, -23.03532101826747, 10.75239953710826), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (8.000015838288707, 3.76649857059908, 18.751976099155204, 10.279138748598198), 0.001 * ra_arcsecond)
    compare(dec.degrees, (20.584000539289498, 19.87762515818775, -23.034049570458144, 10.68696144441038), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (46.72403357148921, 2.2209469369832675, -5.486505415022833, -2.738407691502772), 0.001 * arcsecond)
    compare(az.degrees, (258.5550717845947, 293.95636637272145, 115.32008451470388, 286.09632001391725), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (46.73947196634785, 2.486840978779398, -5.486505415022833, -2.738407691502772), 0.001 * arcsecond)
    compare(az.degrees, (258.5550717845947, 293.95636637272145, 115.32008451470388, 286.09632001391725), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (46.73962875307822, 2.489379891081043, -5.486505415022833, -2.738407691502772), 0.001 * arcsecond)
    compare(az.degrees, (258.5550717845947, 293.95636637272145, 115.32008451470388, 286.09632001391725), 0.001 * arcsecond)

def test_moon_topocentric_date0():
    jd = JulianDate(tt=2440423.345833333)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.moon).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 12.489955349304847, 0.001 * ra_arcsecond)
    compare(dec.degrees, -5.189705732228467, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 12.463855411284248, 0.001 * ra_arcsecond)
    compare(dec.degrees, -5.022075882873389, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 41.920401350253684, 0.001 * arcsecond)
    compare(az.degrees, 151.19707488767648, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 41.938650930938635, 0.001 * arcsecond)
    compare(az.degrees, 151.19707488767648, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 41.93883624837601, 0.001 * arcsecond)
    compare(az.degrees, 151.19707488767648, 0.001 * arcsecond)

def test_moon_topocentric_date1():
    jd = JulianDate(tt=2448031.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.moon).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 23.663473338211592, 0.001 * ra_arcsecond)
    compare(dec.degrees, 1.2271612889121684, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 23.655459675858097, 0.001 * ra_arcsecond)
    compare(dec.degrees, 1.1749464194370667, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, -47.74510120858724, 0.001 * arcsecond)
    compare(az.degrees, 338.1329529181222, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, -47.74510120858724, 0.001 * arcsecond)
    compare(az.degrees, 338.1329529181222, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, -47.74510120858724, 0.001 * arcsecond)
    compare(az.degrees, 338.1329529181222, 0.001 * arcsecond)

def test_moon_topocentric_date2():
    jd = JulianDate(tt=2451545.0)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.moon).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 14.845679251157012, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.590214641231878, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 14.845444624832783, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.587991888462232, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 36.381265580736006, 0.001 * arcsecond)
    compare(az.degrees, 156.2971102404722, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 36.40348032108538, 0.001 * arcsecond)
    compare(az.degrees, 156.2971102404722, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 36.403705864717196, 0.001 * arcsecond)
    compare(az.degrees, 156.2971102404722, 0.001 * arcsecond)

def test_moon_topocentric_date3():
    jd = JulianDate(tt=2456164.5)
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.moon).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, 16.380804513901403, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.790484629243576, 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, 16.393647715389655, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.818976417687217, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, 28.439387966372543, 0.001 * arcsecond)
    compare(az.degrees, 191.294974272018, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, 28.46951344291743, 0.001 * arcsecond)
    compare(az.degrees, 191.294974272018, 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, 28.46981916998486, 0.001 * arcsecond)
    compare(az.degrees, 191.294974272018, 0.001 * arcsecond)

def test_moon_topocentric_date4():
    jd = JulianDate(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5])
    usno = de405.earth.topos('38.9215 N', '77.0669 W', elevation_m=92.0)

    apparent = usno(jd).observe(de405.moon).apparent()
    ra, dec, distance = apparent.radec()
    compare(ra.hours, (12.489955349304847, 23.663473338211592, 14.845679251157012, 16.380804513901403), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-5.189705732228467, 1.2271612889121684, -11.590214641231878, -21.790484629243576), 0.001 * arcsecond)

    ra, dec, distance = apparent.radec(epoch='date')
    compare(ra.hours, (12.463855411284248, 23.655459675858097, 14.845444624832783, 16.393647715389655), 0.001 * ra_arcsecond)
    compare(dec.degrees, (-5.022075882873389, 1.1749464194370667, -11.587991888462232, -21.818976417687217), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz()
    compare(alt.degrees, (41.920401350253684, -47.74510120858724, 36.381265580736006, 28.439387966372543), 0.001 * arcsecond)
    compare(az.degrees, (151.19707488767648, 338.1329529181222, 156.2971102404722, 191.294974272018), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz('standard')
    compare(alt.degrees, (41.938650930938635, -47.74510120858724, 36.40348032108538, 28.46951344291743), 0.001 * arcsecond)
    compare(az.degrees, (151.19707488767648, 338.1329529181222, 156.2971102404722, 191.294974272018), 0.001 * arcsecond)

    alt, az, distance = apparent.altaz(10.0, 1010.0)
    compare(alt.degrees, (41.93883624837601, -47.74510120858724, 36.403705864717196, 28.46981916998486), 0.001 * arcsecond)
    compare(az.degrees, (151.19707488767648, 338.1329529181222, 156.2971102404722, 191.294974272018), 0.001 * arcsecond)

def test_hipparcos_conversion0():
    line = 'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    star = hipparcos.parse(line)
    compare(star.ra.hours, 2.530301023497941, 0.001 * ra_arcsecond)
    compare(star.dec.degrees, 89.26410950742938, 0.001 * arcsecond)
    ra, dec, distance = de405.earth(tt=2440423.345833333).observe(star).radec()
    compare(ra.hours, 2.528369700052897, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26420852419295, 0.001 * arcsecond)

def test_hipparcos_conversion1():
    line = 'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    star = hipparcos.parse(line)
    compare(star.ra.hours, 2.530301023497941, 0.001 * ra_arcsecond)
    compare(star.dec.degrees, 89.26410950742938, 0.001 * arcsecond)
    ra, dec, distance = de405.earth(tt=2448031.5).observe(star).radec()
    compare(ra.hours, 2.5296910104479484, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26413900274704, 0.001 * arcsecond)

def test_hipparcos_conversion2():
    line = 'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    star = hipparcos.parse(line)
    compare(star.ra.hours, 2.530301023497941, 0.001 * ra_arcsecond)
    compare(star.dec.degrees, 89.26410950742938, 0.001 * arcsecond)
    ra, dec, distance = de405.earth(tt=2451545.0).observe(star).radec()
    compare(ra.hours, 2.5302921836971946, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.26411033462212, 0.001 * arcsecond)

def test_hipparcos_conversion3():
    line = 'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    star = hipparcos.parse(line)
    compare(star.ra.hours, 2.530301023497941, 0.001 * ra_arcsecond)
    compare(star.dec.degrees, 89.26410950742938, 0.001 * arcsecond)
    ra, dec, distance = de405.earth(tt=2456164.5).observe(star).radec()
    compare(ra.hours, 2.53111707532574, 0.001 * ra_arcsecond)
    compare(dec.degrees, 89.2640691384828, 0.001 * arcsecond)

def test_hipparcos_conversion4():
    line = 'H|       11767| |02 31 47.08|+89 15 50.9| 1.97|1|H|037.94614689|+89.26413805| |   7.56|   44.22|  -11.74|  0.39|  0.45|  0.48|  0.47|  0.55|-0.16| 0.05| 0.27|-0.01| 0.08| 0.05| 0.04|-0.12|-0.09|-0.36|  1| 1.22| 11767| 2.756|0.003| 2.067|0.003| | 0.636|0.003|T|0.70|0.00|L| | 2.1077|0.0021|0.014|102| | 2.09| 2.13|   3.97|P|1|A|02319+8915|I| 1| 1| | | |  |   |       |     |     |    |S| |P|  8890|B+88    8 |          |          |0.68|F7:Ib-IIv SB|G\n'
    star = hipparcos.parse(line)
    compare(star.ra.hours, 2.530301023497941, 0.001 * ra_arcsecond)
    compare(star.dec.degrees, 89.26410950742938, 0.001 * arcsecond)
    ra, dec, distance = de405.earth(tt=[2440423.345833333, 2448031.5, 2451545.0, 2456164.5]).observe(star).radec()
    compare(ra.hours, (2.528369700052897, 2.5296910104479484, 2.5302921836971946, 2.53111707532574), 0.001 * ra_arcsecond)
    compare(dec.degrees, (89.26420852419295, 89.26413900274704, 89.26411033462212, 89.2640691384828), 0.001 * arcsecond)

