# -*- coding: utf-8 -*-

from numpy import array
from skyfield import api
from skyfield.api import EarthSatellite, load
from skyfield.constants import AU_KM, AU_M
from skyfield.sgp4lib import TEME_to_ITRF
from skyfield.timelib import julian_date

line1 = '1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993'
line2 = '2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106'

# Here are numbers from HORIZONS, which I copied into the test below:
#
#Ephemeris / WWW_USER Wed Jul  4 19:16:45 2018 Pasadena, USA      / Horizons
#...
#2458303.500000000 = A.D. 2018-Jul-04 00:00:00.0000 TDB
# X = 2.633404251158200E-05 Y = 1.015087620439817E-05 Z = 3.544778677556393E-05
# VX=-1.751248694205384E-03 VY= 4.065407557020968E-03 VZ= 1.363540232307603E-04
#2458304.500000000 = A.D. 2018-Jul-05 00:00:00.0000 TDB
# X =-2.136440257814821E-05 Y =-2.084170814514480E-05 Z =-3.415494123796893E-05
# VX= 2.143876266215405E-03 VY=-3.752167957502106E-03 VZ= 9.484159290242074E-04

# TODO: try with array of dates

def test_iss_against_horizons():
    ts = api.load.timescale()
    s = EarthSatellite(line1, line2)

    hp = array([
        [2.633404251158200E-5, 1.015087620439817E-5, 3.544778677556393E-5],
        [-2.136440257814821E-5, -2.084170814514480E-5, -3.415494123796893E-5],
    ]).T
    hv = array([
        [-1.751248694205384E-3, 4.065407557020968E-3, 1.363540232307603E-4],
        [2.143876266215405E-3, -3.752167957502106E-3, 9.484159290242074E-4],
    ]).T

    two_meters = 2.0 / AU_M
    three_km_per_hour = 3.0 * 24.0 / AU_KM

    t = ts.tdb(2018, 7, 4)
    p = s.at(t)
    assert abs(p.position.au - hp[:,0]).max() < two_meters
    assert abs(p.velocity.au_per_d - hv[:,0]).max() < three_km_per_hour

    t = ts.tdb(2018, 7, [4, 5])
    p = s.at(t)
    assert abs(p.position.au - hp).max() < two_meters
    assert abs(p.velocity.au_per_d - hv).max() < three_km_per_hour

# The following tests are based on the text of
# http://www.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf

appendix_c_example = """\
TEME EXAMPLE
1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667
"""

from ..constants import DEG2RAD

arcminute = DEG2RAD / 60.0
arcsecond = arcminute / 60.0
seconds_per_day = 86400.0

# Note that the following test is based specifically on Revision 2 of
# "Revisiting Spacetrack Report #3" AIAA 2006-6753 (earlier versions of
# the PDF use different numbers):
#
# http://ww.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf

def test_appendix_c_conversion_from_TEME_to_ITRF():
    rTEME = array([5094.18016210, 6127.64465950, 6380.34453270])
    vTEME = array([-4.746131487, 0.785818041, 5.531931288])
    vTEME = vTEME * 24.0 * 60.0 * 60.0  # km/s to km/day

    jd_utc = julian_date(2004, 4, 6, 7, 51, 28.386)
    d_ut1 = -0.439961
    jd_ut1 = jd_utc + d_ut1 / 86400.0

    xp = -0.140682 * arcsecond
    yp = 0.333309 * arcsecond

    rITRF, vITRF = TEME_to_ITRF(jd_ut1, rTEME, vTEME, xp, yp)

    epsilon = 5e-8  # Why not 1e-8, which would match all of their digits?

    assert abs(-1033.47938300 - rITRF[0]) < epsilon
    assert abs(+7901.29527540 - rITRF[1]) < epsilon
    assert abs(+6380.35659580 - rITRF[2]) < epsilon

    vITRF_per_second = vITRF / seconds_per_day

    epsilon = 7e-8  # Why not 1e-9, which would match all of their digits?

    assert abs(-3.225636520 - vITRF_per_second[0]) < epsilon
    assert abs(-2.872451450 - vITRF_per_second[1]) < epsilon
    assert abs(+5.531924446 - vITRF_per_second[2]) < epsilon

def test_appendix_c_satellite():
    lines = appendix_c_example.splitlines()
    ts = api.load.timescale()
    sat = EarthSatellite(lines[1], lines[2], lines[0], ts)
    t = ts.tt_jd(sat.epoch.whole + 3.0, sat.epoch.tt_fraction)

    # First, a crucial sanity check (which is, technically, a test of
    # the `sgp4` package and not of Skyfield): are the right coordinates
    # being produced by our Python SGP4 propagator for this satellite?

    rTEME, vTEME, error = sat._position_and_velocity_TEME_km(t)

    # TODO: This used to be accurate to within 1e-8 but lost precision
    # with the move to SGP4 2.0.  Is the difference an underlying change
    # in the algorithm and its results?  Or something else?
    epsilon = 1e-4
    assert abs(-9060.47373569 - rTEME[0]) < epsilon
    assert abs(4658.70952502 - rTEME[1]) < epsilon
    assert abs(813.68673153 - rTEME[2]) < epsilon

    # TODO: Similar to the above, this used to be 1e-9.
    epsilon = 1e-8
    assert abs(-2.232832783 - vTEME[0]) < epsilon
    assert abs(-4.110453490 - vTEME[1]) < epsilon
    assert abs(-3.157345433 - vTEME[2]) < epsilon

def test_epoch_date():
    # Example from https://celestrak.com/columns/v04n03/
    s = appendix_c_example.replace('00179.78495062', '98001.00000000')
    lines = s.splitlines()
    sat = EarthSatellite(lines[1], lines[2], lines[0])
    assert sat.epoch.utc_jpl() == 'A.D. 1998-Jan-01 00:00:00.0000 UTC'

def test_target_number():
    s = EarthSatellite(line1, line2)
    assert s.target == -125544

def test_is_sunlit():
    # Yes, a positionlib method; but it made sense to test it here.
    ts = api.load.timescale()
    t = ts.utc(2018, 7, 3, 0, range(0, 60, 10))
    s = EarthSatellite(line1, line2)
    eph = load('de421.bsp')
    expected = [True, False, False, False, True, True]
    assert list(s.at(t).is_sunlit(eph)) == expected

    # What if we observe from a topos rather than the geocenter?
    topos = api.Topos('40.8939 N', '83.8917 W')
    assert list((s - topos).at(t).is_sunlit(eph)) == expected

def test_is_venus_behind_earth():
    # Like the previous test: a satellite-focused positionlib method.
    # Just for fun, we ask whether the Sun is behind the earth, so this
    # measures the same celestial circumstance as the previous test.
    ts = api.load.timescale()
    t = ts.utc(2018, 7, 3, 0, range(0, 60, 10))
    s = EarthSatellite(line1, line2)
    eph = load('de421.bsp')
    expected = [False, True, True, True, False, False]
    p = (eph['earth'] + s).at(t).observe(eph['sun']).apparent()
    assert list(p.is_behind_earth()) == expected

def test_is_another_satellite_behind_earth():
    # See if the method works with a pure geometric difference.
    ts = api.load.timescale()
    t = ts.utc(2018, 7, 3, 0, range(0, 60, 10))
    s = EarthSatellite(line1, line2)
    # The "other satellite" is fictitious: the ISS offset by one day.
    s2 = EarthSatellite(line1.replace('184.80969102', '185.80969102'), line2)
    expected = [True, True, True, True, True, True]
    p = (s - s2).at(t)
    assert list(p.is_behind_earth()) == expected
