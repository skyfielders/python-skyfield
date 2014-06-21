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

def test_mercury_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 1.3278115470600746, 0.5 * meter)

    a = e.observe(de405.mercury)
    ra, dec, distance = a.radec()
    compare(ra.hours, 7.905384000977572, 0.001 * ra_arcsecond)
    compare(dec.degrees, 22.332364359841474, 0.001 * arcsecond)

def test_venus_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 0.9646045654448725, 0.5 * meter)

    a = e.observe(de405.venus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 4.966946050917652, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.210417323471006, 0.001 * arcsecond)

def test_mars_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 0.5912188976380217, 0.5 * meter)

    a = e.observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 16.0296606272219, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.12731030858147, 0.001 * arcsecond)

def test_jupiter_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 5.841600319231745, 0.5 * meter)

    a = e.observe(de405.jupiter)
    ra, dec, distance = a.radec()
    compare(ra.hours, 12.104091505864655, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.6513409058207986, 0.001 * arcsecond)

def test_saturn_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 9.382032444401027, 0.5 * meter)

    a = e.observe(de405.saturn)
    ra, dec, distance = a.radec()
    compare(ra.hours, 2.462774885242021, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.045819985925936, 0.001 * arcsecond)

def test_uranus_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 18.75197906203834, 0.5 * meter)

    a = e.observe(de405.uranus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 12.087167068351333, 0.001 * ra_arcsecond)
    compare(dec.degrees, 0.20723926118363262, 0.001 * arcsecond)

def test_neptune_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 29.83221264621946, 0.5 * meter)

    a = e.observe(de405.neptune)
    ra, dec, distance = a.radec()
    compare(ra.hours, 15.637210587139663, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.67999613660563, 0.001 * arcsecond)

def test_pluto_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 32.312971776632494, 0.5 * meter)

    a = e.observe(de405.pluto)
    ra, dec, distance = a.radec()
    compare(ra.hours, 12.015311208821212, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.620557180992584, 0.001 * arcsecond)

def test_sun_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 1.0160878650466754, 0.5 * meter)

    a = e.observe(de405.sun)
    ra, dec, distance = a.radec()
    compare(ra.hours, 8.03008088792976, 0.001 * ra_arcsecond)
    compare(dec.degrees, 20.496475643233936, 0.001 * arcsecond)

def test_moon_astrometric_0():
    jd = JulianDate(tt=2440423.345833333)

    e = de405.earth(jd)
    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.0026034424248855136, 0.5 * meter)

    a = e.observe(de405.moon)
    ra, dec, distance = a.radec()
    compare(ra.hours, 12.472463241145173, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.546618838171283, 0.001 * arcsecond)

def test_mercury_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 0.6507044512046538, 0.5 * meter)

    a = e.observe(de405.mercury)
    ra, dec, distance = a.radec()
    compare(ra.hours, 2.4704717994133576, 0.001 * ra_arcsecond)
    compare(dec.degrees, 11.250132844930501, 0.001 * arcsecond)

def test_venus_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 1.0711674186789975, 0.5 * meter)

    a = e.observe(de405.venus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 1.161811406279447, 0.001 * ra_arcsecond)
    compare(dec.degrees, 5.32829157368082, 0.001 * arcsecond)

def test_mars_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 1.4302506796029135, 0.5 * meter)

    a = e.observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 23.545034875459514, 0.001 * ra_arcsecond)
    compare(dec.degrees, -4.882249043221036, 0.001 * arcsecond)

def test_jupiter_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 5.913287883102949, 0.5 * meter)

    a = e.observe(de405.jupiter)
    ra, dec, distance = a.radec()
    compare(ra.hours, 6.765154678701348, 0.001 * ra_arcsecond)
    compare(dec.degrees, 23.17039770012201, 0.001 * arcsecond)

def test_saturn_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 9.420484451056101, 0.5 * meter)

    a = e.observe(de405.saturn)
    ra, dec, distance = a.radec()
    compare(ra.hours, 19.814248756112033, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.93339019805076, 0.001 * arcsecond)

def test_uranus_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 18.62241700929518, 0.5 * meter)

    a = e.observe(de405.uranus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 18.668551452013403, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.437331340689155, 0.001 * arcsecond)

def test_neptune_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 29.490001740438892, 0.5 * meter)

    a = e.observe(de405.neptune)
    ra, dec, distance = a.radec()
    compare(ra.hours, 19.03623522579387, 0.001 * ra_arcsecond)
    compare(dec.degrees, -21.792864018500975, 0.001 * arcsecond)

def test_pluto_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 28.707485955458118, 0.5 * meter)

    a = e.observe(de405.pluto)
    ra, dec, distance = a.radec()
    compare(ra.hours, 15.216302246424346, 0.001 * ra_arcsecond)
    compare(dec.degrees, -1.3346560528819575, 0.001 * arcsecond)

def test_sun_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 1.0118605934887042, 0.5 * meter)

    a = e.observe(de405.sun)
    ra, dec, distance = a.radec()
    compare(ra.hours, 3.776110727862678, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.907832379364578, 0.001 * arcsecond)

def test_moon_astrometric_1():
    jd = JulianDate(tt=2448031.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.002481509229659986, 0.5 * meter)

    a = e.observe(de405.moon)
    ra, dec, distance = a.radec()
    compare(ra.hours, 23.67644381740948, 0.001 * ra_arcsecond)
    compare(dec.degrees, 1.8587554901327863, 0.001 * arcsecond)

def test_mercury_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 1.4155249674526948, 0.5 * meter)

    a = e.observe(de405.mercury)
    ra, dec, distance = a.radec()
    compare(ra.hours, 18.13892977357885, 0.001 * ra_arcsecond)
    compare(dec.degrees, -24.420324941080732, 0.001 * arcsecond)

def test_venus_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 1.1376890757925104, 0.5 * meter)

    a = e.observe(de405.venus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 15.993350650200568, 0.001 * ra_arcsecond)
    compare(dec.degrees, -18.451653207795232, 0.001 * arcsecond)

def test_mars_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 1.8496039270835372, 0.5 * meter)

    a = e.observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 22.034936616343344, 0.001 * ra_arcsecond)
    compare(dec.degrees, -13.180707411034978, 0.001 * arcsecond)

def test_jupiter_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 4.621126565890217, 0.5 * meter)

    a = e.observe(de405.jupiter)
    ra, dec, distance = a.radec()
    compare(ra.hours, 1.5913207023268698, 0.001 * ra_arcsecond)
    compare(dec.degrees, 8.595887646396902, 0.001 * arcsecond)

def test_saturn_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 8.652750126001484, 0.5 * meter)

    a = e.observe(de405.saturn)
    ra, dec, distance = a.radec()
    compare(ra.hours, 2.584400980536592, 0.001 * ra_arcsecond)
    compare(dec.degrees, 12.616288735770388, 0.001 * arcsecond)

def test_uranus_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 20.727159134679393, 0.5 * meter)

    a = e.observe(de405.uranus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 21.16558686754142, 0.001 * ra_arcsecond)
    compare(dec.degrees, -17.018831731314233, 0.001 * arcsecond)

def test_neptune_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 31.02449192035449, 0.5 * meter)

    a = e.observe(de405.neptune)
    ra, dec, distance = a.radec()
    compare(ra.hours, 20.362841834121518, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.212425239376326, 0.001 * arcsecond)

def test_pluto_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 31.064412196006607, 0.5 * meter)

    a = e.observe(de405.pluto)
    ra, dec, distance = a.radec()
    compare(ra.hours, 16.761873062250743, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.39643313463007, 0.001 * arcsecond)

def test_sun_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 0.9833276788862821, 0.5 * meter)

    a = e.observe(de405.sun)
    ra, dec, distance = a.radec()
    compare(ra.hours, 18.752544254682526, 0.001 * ra_arcsecond)
    compare(dec.degrees, -23.033309607967187, 0.001 * arcsecond)

def test_moon_astrometric_2():
    jd = JulianDate(tt=2451545.0)

    e = de405.earth(jd)
    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.0026902029885132865, 0.5 * meter)

    a = e.observe(de405.moon)
    ra, dec, distance = a.radec()
    compare(ra.hours, 14.830020573942209, 0.001 * ra_arcsecond)
    compare(dec.degrees, -10.900635500944452, 0.001 * arcsecond)

def test_mercury_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.mercury(jd)).position.AU)
    compare(distance, 1.1264323486728112, 0.5 * meter)

    a = e.observe(de405.mercury)
    ra, dec, distance = a.radec()
    compare(ra.hours, 9.295934662566733, 0.001 * ra_arcsecond)
    compare(dec.degrees, 16.68579742896488, 0.001 * arcsecond)

def test_venus_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.venus(jd)).position.AU)
    compare(distance, 0.7824924286112767, 0.5 * meter)

    a = e.observe(de405.venus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 7.175585125577371, 0.001 * ra_arcsecond)
    compare(dec.degrees, 19.874130272238094, 0.001 * arcsecond)

def test_mars_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.mars(jd)).position.AU)
    compare(distance, 1.766552316866877, 0.5 * meter)

    a = e.observe(de405.mars)
    ra, dec, distance = a.radec()
    compare(ra.hours, 13.894324196598355, 0.001 * ra_arcsecond)
    compare(dec.degrees, -12.122808318928705, 0.001 * arcsecond)

def test_jupiter_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.jupiter(jd)).position.AU)
    compare(distance, 5.129958529243068, 0.5 * meter)

    a = e.observe(de405.jupiter)
    ra, dec, distance = a.radec()
    compare(ra.hours, 4.822841055032964, 0.001 * ra_arcsecond)
    compare(dec.degrees, 21.649994488649472, 0.001 * arcsecond)

def test_saturn_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.saturn(jd)).position.AU)
    compare(distance, 10.326368974662916, 0.5 * meter)

    a = e.observe(de405.saturn)
    ra, dec, distance = a.radec()
    compare(ra.hours, 13.628484577191722, 0.001 * ra_arcsecond)
    compare(dec.degrees, -7.659435207931653, 0.001 * arcsecond)

def test_uranus_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.uranus(jd)).position.AU)
    compare(distance, 19.234768680195387, 0.5 * meter)

    a = e.observe(de405.uranus)
    ra, dec, distance = a.radec()
    compare(ra.hours, 0.48916431485643164, 0.001 * ra_arcsecond)
    compare(dec.degrees, 2.3565095329111823, 0.001 * arcsecond)

def test_neptune_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.neptune(jd)).position.AU)
    compare(distance, 28.98411802971635, 0.5 * meter)

    a = e.observe(de405.neptune)
    ra, dec, distance = a.radec()
    compare(ra.hours, 22.252468120719442, 0.001 * ra_arcsecond)
    compare(dec.degrees, -11.504657215501586, 0.001 * arcsecond)

def test_pluto_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.pluto(jd)).position.AU)
    compare(distance, 31.69909782133192, 0.5 * meter)

    a = e.observe(de405.pluto)
    ra, dec, distance = a.radec()
    compare(ra.hours, 18.488351288595236, 0.001 * ra_arcsecond)
    compare(dec.degrees, -19.552190994888846, 0.001 * arcsecond)

def test_sun_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.sun(jd)).position.AU)
    compare(distance, 1.0107820040799866, 0.5 * meter)

    a = e.observe(de405.sun)
    ra, dec, distance = a.radec()
    compare(ra.hours, 10.268162490439074, 0.001 * ra_arcsecond)
    compare(dec.degrees, 10.751933902906117, 0.001 * arcsecond)

def test_moon_astrometric_3():
    jd = JulianDate(tt=2456164.5)

    e = de405.earth(jd)
    distance = length_of((e - de405.moon(jd)).position.AU)
    compare(distance, 0.0024739078649309238, 0.5 * meter)

    a = e.observe(de405.moon)
    ra, dec, distance = a.radec()
    compare(ra.hours, 16.39102815233177, 0.001 * ra_arcsecond)
    compare(dec.degrees, -20.93676001523414, 0.001 * arcsecond)

