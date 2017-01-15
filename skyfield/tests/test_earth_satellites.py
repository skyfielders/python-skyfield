# -*- coding: utf-8 -*-

import sys
from datetime import datetime, timedelta
from numpy import array
from skyfield import api
from skyfield.sgp4lib import EarthSatellite, TEME_to_ITRF
from skyfield.timelib import utc

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

def ts():
    yield api.load.timescale()

def iss_transit():
    for line in heavens_above_transits.splitlines():
        fields = line.split()
        dt = datetime.strptime('2013 {0} {1} {6}'.format(*fields),
                               '%Y %d %b %H:%M:%S').replace(tzinfo=utc)
        altitude = float(fields[7][:-1])
        yield dt, altitude

def _OFF_test_iss_altitude_computed_with_bcrs(iss_transit):
    # This test has been disabled because I no longer intend to support
    # the complexity of implicit operations between BCRS positions and
    # GCRS positions.
    dt, their_altitude = iss_transit

    cst = timedelta(hours=-6) #, minutes=1)
    dt = dt - cst
    t = api.load.timescale(delta_t=67.2091).utc(dt)

    lines = iss_tle.splitlines()
    s = EarthSatellite(lines, None)
    earth = api.load('de421.bsp')['earth']
    lake_zurich = earth.topos(latitude_degrees=42.2, longitude_degrees=-88.1)

    # Compute using Solar System coordinates:

    alt, az, d = lake_zurich.at(t).observe(s).altaz()
    print(dt, their_altitude, alt.degrees, their_altitude - alt.degrees)
    assert abs(alt.degrees - their_altitude) < 2.5  # TODO: tighten this up?

def test_iss_altitude_computed_with_gcrs(iss_transit):
    dt, their_altitude = iss_transit

    cst = timedelta(hours=-6) #, minutes=1)
    dt = dt - cst
    t = api.load.timescale(delta_t=67.2091).utc(dt)

    lines = iss_tle.splitlines()
    s = EarthSatellite(lines, None)
    lake_zurich = api.Topos(latitude_degrees=42.2, longitude_degrees=-88.1)

    alt, az, d = lake_zurich.at(t).observe(s).altaz()
    print(dt, their_altitude, alt.degrees, their_altitude - alt.degrees)
    assert abs(alt.degrees - their_altitude) < 2.5  # TODO: tighten this up?

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
    sat = EarthSatellite(lines, None)

    ts = api.load.timescale()
    jd_epoch = sat._sgp4_satellite.jdsatepoch
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
    sat = EarthSatellite(s.splitlines(), None)
    assert sat.epoch.utc_jpl() == 'A.D. 1998-Jan-01 00:00:00.0000 UT'
