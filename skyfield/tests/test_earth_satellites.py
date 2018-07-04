# -*- coding: utf-8 -*-

from numpy import array
from skyfield import api
from skyfield.api import EarthSatellite
from skyfield.constants import AU_KM, AU_M
from skyfield.sgp4lib import TEME_to_ITRF

iss_tle0 = """\
1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993
2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106
"""

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
    s = EarthSatellite(*iss_tle0.splitlines())

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
    assert (abs(p.position.au - hp[:,0]) < two_meters).all()
    assert (abs(p.velocity.au_per_d - hv[:,0]) < three_km_per_hour).all()

    t = ts.tdb(2018, 7, [4, 5])
    p = s.at(t)
    assert (abs(p.position.au - hp) < two_meters).all()
    assert (abs(p.velocity.au_per_d - hv) < three_km_per_hour).all()

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
second = 1.0 / (24.0 * 60.0 * 60.0)

# Note that the following test is based specifically on Revision 2 of
# "Revisiting Spacetrack Report #3" AIAA 2006-6753 (earlier versions of
# the PDF use different numbers):
#
# http://ww.celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf

def test_appendix_c_conversion_from_TEME_to_ITRF():
    rTEME = array([5094.18016210, 6127.64465950, 6380.34453270])
    vTEME = array([-4.746131487, 0.785818041, 5.531931288])
    vTEME = vTEME * 24.0 * 60.0 * 60.0  # km/s to km/day

    ts = api.load.timescale()
    jd_ut1 = ts.tt(2004, 4, 6, 7, 51, 28.386 - 0.439961).tt

    xp = -0.140682 * arcsecond
    yp = 0.333309 * arcsecond

    rITRF, vITRF = TEME_to_ITRF(jd_ut1, rTEME, vTEME, xp, yp)

    meter = 1e-3

    assert abs(-1033.47938300 - rITRF[0]) < 0.1 * meter
    assert abs(+7901.29527540 - rITRF[1]) < 0.1 * meter
    assert abs(+6380.35659580 - rITRF[2]) < 0.1 * meter

    vITRF_per_second = vITRF * second

    assert abs(-3.225636520 - vITRF_per_second[0]) < 1e-4 * meter
    assert abs(-2.872451450 - vITRF_per_second[1]) < 1e-4 * meter
    assert abs(+5.531924446 - vITRF_per_second[2]) < 1e-4 * meter

def test_appendix_c_satellite():
    lines = appendix_c_example.splitlines()
    sat = EarthSatellite(lines[1], lines[2], lines[0])

    ts = api.load.timescale()
    jd_epoch = sat.model.jdsatepoch
    three_days_later = jd_epoch + 3.0
    offset = ts.tt(jd=three_days_later)._utc_float() - three_days_later
    t = ts.tt(jd=three_days_later - offset)

    # First, a crucial sanity check (which is, technically, a test of
    # the `sgp4` package and not of Skyfield): are the right coordinates
    # being produced by our Python SGP4 propagator for this satellite?

    rTEME, vTEME, error = sat._position_and_velocity_TEME_km(t)

    assert abs(-9060.47373569 - rTEME[0]) < 1e-8
    assert abs(4658.70952502 - rTEME[1]) < 1e-8
    assert abs(813.68673153 - rTEME[2]) < 1e-8

    assert abs(-2.232832783 - vTEME[0]) < 1e-9
    assert abs(-4.110453490 - vTEME[1]) < 1e-9
    assert abs(-3.157345433 - vTEME[2]) < 1e-9

def test_epoch_date():
    # Example from https://celestrak.com/columns/v04n03/
    s = appendix_c_example.replace('00179.78495062', '98001.00000000')
    lines = s.splitlines()
    sat = EarthSatellite(lines[1], lines[2], lines[0])
    assert sat.epoch.utc_jpl() == 'A.D. 1998-Jan-01 00:00:00.0000 UT'
