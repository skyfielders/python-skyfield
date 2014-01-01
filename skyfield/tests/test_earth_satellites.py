# -*- coding: utf-8 -*-

import pytest
import sys
from datetime import datetime, timedelta
from numpy import array
from skyfield.planets import earth
from skyfield.sgp4lib import EarthSatellite, TEME_to_ITRF
from skyfield.timescales import JulianDate, julian_date

iss_tle = ("""\
ISS (ZARYA)             \n\
1 25544U 98067A   13330.58127943  .00000814  00000-0  21834-4 0  1064\n\
2 25544  51.6484  23.7537 0001246  74.1647  18.7420 15.50540527859894\n\
""")

heavens_above_transits = """\
26 Nov	-1.7	04:55:55	29°	NNE	04:55:55	29°	NNE	04:58:45	10°	E	visible
27 Nov	0.1	04:09:07	12°	ENE	04:09:07	12°	ENE	04:09:25	10°	ENE	visible
27 Nov	-3.4	05:42:00	26°	WNW	05:43:45	86°	SW	05:47:05	10°	SE	visible
28 Nov	-2.6	04:55:15	52°	ENE	04:55:15	52°	ENE	04:58:07	10°	ESE	visible
29 Nov	0.1	04:08:35	13°	E	04:08:35	13°	E	04:08:58	10°	E	visible
29 Nov	-2.2	05:41:28	25°	WSW	05:42:30	31°	SW	05:45:28	10°	SSE	visible
30 Nov	-1.9	04:54:52	33°	SSE	04:54:52	33°	SSE	04:56:55	10°	SE	visible
01 Dec	-0.9	05:41:13	12°	SW	05:41:13	12°	SW	05:42:15	10°	SSW	visible
02 Dec	-0.4	04:54:46	10°	S	04:54:46	10°	S	04:54:49	10°	S	visible
"""
if sys.version_info < (3,):
    heavens_above_transits = heavens_above_transits.decode('utf-8')

@pytest.fixture(params=heavens_above_transits.splitlines())
def iss_transit(request):
    line = request.param
    fields = line.split()
    dt = datetime.strptime('2013 {0} {1} {6}'.format(*fields),
                           '%Y %d %b %H:%M:%S')
    altitude = float(fields[7][:-1])
    return dt, altitude

def test_iss_altitude(iss_transit):
    dt, their_altitude = iss_transit

    cst = timedelta(hours=-6) #, minutes=1)
    dt = dt - cst
    jd = JulianDate(dt)

    lines = iss_tle.splitlines()
    s = EarthSatellite(lines, earth)
    lake_zurich = earth.topos('88.1 W', '42.2 N')
    alt, az, d = lake_zurich(jd).observe(s).altaz()
    print(dt, their_altitude, alt.degrees(), their_altitude - alt.degrees())
    assert abs(alt.degrees() - their_altitude) < 2.5  # TODO: tighten this up?

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

def test_appendix_c_conversion_from_TEME_to_ITRF():
    rTEME = array([5094.18016210, 6127.64465950, 6380.34453270])
    vTEME = array([-4.746131487, 0.785818041, 5.531931288])
    vTEME = vTEME * 24.0 * 60.0 * 60.0  # km/s to km/day

    jd_ut1 = julian_date(2004, 4, 6, 7, 51, 28.386 - 0.439961)
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
    sat = EarthSatellite(lines, earth)

    jd_epoch = sat._sgp4_satellite.jdsatepoch
    three_days_later = jd_epoch + 3.0
    jd = JulianDate(ut1=three_days_later)

    # First, a crucial sanity check (which is, technically, a test of
    # the `sgp4` package and not of Skyfield): are the right coordinates
    # being produced by our Python SGP4 propagator for this satellite?

    rTEME, vTEME = sat._position_and_velocity_TEME_km(jd)

    assert abs(-9060.47373569 - rTEME[0]) < 1e-8
    assert abs(4658.70952502 - rTEME[1]) < 1e-8
    assert abs(813.68673153 - rTEME[2]) < 1e-8

    assert abs(-2.232832783 - vTEME[0]) < 1e-9
    assert abs(-4.110453490 - vTEME[1]) < 1e-9
    assert abs(-3.157345433 - vTEME[2]) < 1e-9
