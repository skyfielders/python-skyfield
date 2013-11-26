from numpy import array, round
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from skyfield.constants import KM_AU
from skyfield.planets import earth, mars
from skyfield.timescales import JulianDate, julian_date

iss_tle = (
    'ISS (ZARYA)             \n',
    '1 25544U 98067A   13329.81056895  .00017555  00000-0  30885-3 0  1015\n',
    '2 25544  51.6512  27.5741 0001256 103.0346   4.9903 15.50551370859772\n',
    )

def test_iss_altitude():
    #
    # Transit times and altitudes predicted by spaceweather.com.
    #
    # Lake Zurich, Illin (42.2 lat., -88.1 long.)
    # Nov 24
    # ISS 06:30:23 am WNW 06:33:43am 63 deg -3.7 (very bright)
    #
    cst = -6
    line0, line1, line2 = iss_tle
    satellite = twoline2rv(line1, line2, wgs72)
    iss_position, iss_velocity = satellite.propagate(
        2013, 11, 24, 6 - cst, 33., 43.)
    iss_position = array(iss_position) * KM_AU
    jd = JulianDate(ut1=julian_date(
        2013, 11, 24, (6 - cst) + 33/60. + 43/3600.))
    lake_zurich = earth.topos('88.1 W', '42.2 N')
    earthpos = earth(jd)
    lake = (lake_zurich(jd).position - earthpos.position)
    diff = lake_zurich(jd).observe(mars)
    diff.position = (iss_position) - lake
    diff.lighttime = 0.0000000001
    diff.observer.velocity = array([0.000000001,0.,0.])
    alt, az, d = diff.apparent().altaz()
    assert round(alt.degrees()) == 63.
