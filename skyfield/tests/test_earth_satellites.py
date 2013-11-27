# -*- coding: utf-8 -*-

import pytest
from datetime import datetime, timedelta
from numpy import array
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from skyfield.constants import KM_AU
from skyfield.planets import earth, mars
from skyfield.timescales import JulianDate

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
""".decode('utf-8')

@pytest.fixture(params=heavens_above_transits.splitlines())
def iss_transit(request):
    line = request.param
    fields = line.split()
    dt = datetime.strptime('2013 {0} {1} {6}'.format(*fields),
                           '%Y %d %b %H:%M:%S')
    altitude = float(fields[7][:-1])
    return dt, altitude

def test_iss_altitude(iss_transit):
    dt, altitude = iss_transit
    cst = timedelta(hours=-6) #, minutes=1)
    dt = dt - cst
    line0, line1, line2 = iss_tle.splitlines()
    satellite = twoline2rv(line1, line2, wgs72)
    iss_position, iss_velocity = satellite.propagate(
        dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    iss_position = array(iss_position) * KM_AU
    jd = JulianDate(dt)
    lake_zurich = earth.topos('88.1 W', '42.2 N')
    earthpos = earth(jd)
    lake = (lake_zurich(jd).position - earthpos.position)
    diff = lake_zurich(jd).observe(mars)
    diff.position = (iss_position) - lake
    diff.lighttime = 0.0000000000000001
    diff.observer.velocity = array([0.000000000000001,0.,0.])
    alt, az, d = diff.apparent().altaz()
    print dt, altitude, alt.degrees()
    assert abs(alt.degrees() - altitude) < 2.0  # TODO: tighten this up?
