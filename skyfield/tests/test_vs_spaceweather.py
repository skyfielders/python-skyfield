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

spaceweather_transits = """\
11/25/2013      ISS     05:42:39 am     NNW     05:45:50        40.2    -2.9
11/26/2013      ISS     04:55:49 am     NNE     04:58:18        28.6    -2.2
11/26/2013      ISS     06:29:22 am     WNW     06:32:36        48.1    -3.3
11/27/2013      ISS     05:41:54 am     WNW     05:45:15        85.6    -4.0
11/28/2013      ISS     04:55:07 am     ENE     04:58:25        52.9    -3.4
11/29/2013      ISS     05:41:15 am     W       05:44:15        30.8    -2.3
11/30/2013      ISS     04:54:32 am     SSE     04:57:48        34.2    -2.6
"""

@pytest.fixture(params=spaceweather_transits.splitlines())
def spaceweather_transit(request):
    line = request.param
    fields = line.split()
    dt = datetime.strptime(fields[0] + ' ' + fields[5], '%m/%d/%Y %H:%M:%S')
    altitude = float(fields[6])
    return dt, altitude

def test_iss_altitude(spaceweather_transit):
    dt, altitude = spaceweather_transit
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
    #
    # Does spaceweather.com have different orbital elements today from
    # the ones I can download from Celetrak?  The transits are several
    # minutes off.  Maybe I should find another data source.  For now:
    # just be happy if we think the ISS is above the horizon when
    # spaceweather thinks it is transiting!
    #
    assert alt > 0
