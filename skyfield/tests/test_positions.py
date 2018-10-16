from __future__ import print_function

import numpy as np
from skyfield import api
from skyfield.earthlib import earth_rotation_angle
from skyfield.positionlib import ICRF
from skyfield.starlib import Star
from skyfield.api import Topos, utc
from datetime import datetime

def test_J2000_ecliptic_coordinates_with_and_without_a_time_array():
    p0 = ICRF((1,0,0))
    p1 = ICRF((0,1,0))
    p2 = ICRF(((1, 0),
               (0, 1),
               (0, 0)))

    lat0, lon0, distance0 = p0.ecliptic_latlon(epoch=None)
    lat1, lon1, distance1 = p1.ecliptic_latlon(epoch=None)
    lat2, lon2, distance2 = p2.ecliptic_latlon(epoch=None)

    assert lat2.degrees[0] == lat0.degrees
    assert lat2.degrees[1] == lat1.degrees

    assert lon2.degrees[0] == lon0.degrees
    assert lon2.degrees[1] == lon1.degrees

    assert distance2.au[0] == distance0.au
    assert distance2.au[1] == distance1.au

def test_dynamic_ecliptic_coordinates_with_and_without_a_time_array():
    ts = api.load.timescale()
    t = ts.utc(1980)

    p0 = ICRF((1,0,0))
    p1 = ICRF((0,1,0))
    p2 = ICRF(((1, 0),
               (0, 1),
               (0, 0)))

    lat0, lon0, distance0 = p0.ecliptic_latlon(epoch=t)
    lat1, lon1, distance1 = p1.ecliptic_latlon(epoch=t)
    lat2, lon2, distance2 = p2.ecliptic_latlon(epoch=t)

    assert lat2.degrees[0] == lat0.degrees
    assert lat2.degrees[1] == lat1.degrees

    assert lon2.degrees[0] == lon0.degrees
    assert lon2.degrees[1] == lon1.degrees

    assert distance2.au[0] == distance0.au
    assert distance2.au[1] == distance1.au

# Test that the CIRS coordinate of the TIO is consistent with the Earth Rotation Angle
# This is mostly an internal consistency check
def test_cirs_era():

    ts = api.load.timescale()
    st = ts.utc(year=np.arange(1951, 2051))

    planets = api.load('de421.bsp')
    pos = planets['earth'] + api.Topos(longitude_degrees=0.0, latitude_degrees=0.0)

    # Get the TIO
    tio = pos.at(st).from_altaz(alt_degrees=90, az_degrees=180)

    # Get the TIOs RA in CIRS coordinates, and the Earth Rotation Angle
    tio_ra, tio_dec, _ = tio.cirs_radec(st)
    era = 360.0 * earth_rotation_angle(st.ut1)

    tol = (1e-8 / 3600.0)  # 10 nano arc-second precision
    assert np.allclose(tio_ra._degrees, era, rtol=0.0, atol=tol)
    assert np.allclose(tio_dec._degrees, 0.0, rtol=0.0, atol=tol)


# Check a line of points along the terrestrial prime meridian all have the same
# CIRS RA, and that their declinations are correct.
def test_cirs_meridian():

    ts = api.load.timescale()
    st = ts.utc(year=2051)

    planets = api.load('de421.bsp')
    pos = planets['earth'] + api.Topos(longitude_degrees=0.0, latitude_degrees=0.0)

    # Get a series of points along the meridian
    alt = np.arange(1, 90)
    meridian = pos.at(st).from_altaz(alt_degrees=alt, az_degrees=0.0)

    # Get the TIOs RA in CIRS coordinates, and the Earth Rotation Angle
    md_ra, md_dec, _ = meridian.cirs_radec(st)
    era = 360.0 * earth_rotation_angle(st.ut1)

    tol = (1e-7 / 3600.0)  # 100 nano arc-second precision
    assert np.allclose(md_ra._degrees, era, rtol=0.0, atol=tol)
    assert np.allclose(md_dec._degrees, 90 - alt, rtol=0.0, atol=tol)


# Check a set of positions and times against results calculated externally
# using the IAU SOFA library (20180130 release). For reference the code used
# was:
#
# #include <stdio.h>
# #include <sofa.h>
# #include <math.h>
#
# int main(int argc, char ** argv) {
#
#     // Test data as RA, DEC, TDB. Positions ICRS(deg), time in JD.
#     double test_data[3][3] = {
#         {45.0, 46.0, 2458327},
#         {200.0, -22.0, 2458327},
#         {45.0, 46.0, 2459327}
#     };
#
#     for(int i = 0; i < 3; i++) {
#         double ra_icrs = test_data[i][0] / 180.0 * M_PI;
#         double dec_icrs = test_data[i][1] / 180.0 * M_PI;
#         double jd_tdb = test_data[i][2];
#
#         double ra_cirs, dec_cirs;
#         double eo;
#
#         iauAtci13(ra_icrs, dec_icrs, 0.0, 0.0, 0.0, 0.0, jd_tdb, 0.0,
#                   &ra_cirs, &dec_cirs, &eo);
#
#         printf("%.12f  %.12f\n", ra_cirs / (2 * M_PI) * 360,
#                dec_cirs / (2 * M_PI) * 360);
#     }
# }
def test_cirs_sofa():
    ts = api.load.timescale()
    earth = api.load('de421.bsp')['earth']

    test_data = [
        [45.0, 46.0, 2458327],
        [200.0, -22.0, 2458327],
        [45.0, 46.0, 2459327]
    ]

    # Results output by SOFA. Calculated using the source code above.
    sofa_results = [
        [45.074343838325, 46.067831092355],
        [200.013551320030,  -22.096008994214],
        [45.077698288877,  46.082296559677]
    ]

    tol = 1e-5 / 3600.0  # 10 micro arc-seconds

    for ((ra_icrs, dec_icrs, tdb), (ra_sofa, dec_sofa)) in zip(test_data, sofa_results):
        ss = Star(ra_hours=(ra_icrs / 15.0), dec_degrees=dec_icrs)
        st = ts.tdb(jd=tdb)
        ra_cirs, dec_cirs, _ = earth.at(st).observe(ss).apparent().cirs_radec(st)

        assert np.allclose(ra_cirs._degrees, ra_sofa, rtol=0.0, atol=tol)
        assert np.allclose(dec_cirs._degrees, dec_sofa, rtol=0.0, atol=tol)

def test_east_north_up():
    ts = api.load.timescale()

    # create test reference points on the earth as lat/lon, degrees
    test_sites = [ [ 0.0, 0.0, 0.0 ], # lat/lon,elevation at equator, prime meridian
                   [ 0.0, 0.0, 0.0 ],
                   [ 0.0, 0.0, 0.0 ],
                   [ 0.0, 0.0, 0.0 ],
                  ]
    

    # create a list of test targets to view from each test site
    test_targets = [ [ 0.0, 0.0, 100e6], # lat,lon,elevation
                     [ 0.0, 0.0, 100.0],
                     [ 90.0, 0.0, 10e12],
                     [ 90.0, 0.0, 0.0 ],
                     ]

    # truth results expected in east north, up, meters
    from skyfield.constants import ERAD, IERS_2010_INVERSE_EARTH_FLATTENING
    m_rad = ERAD # meridional radius
    p_rad = ERAD * ( 1-1.0/IERS_2010_INVERSE_EARTH_FLATTENING) # polar radius
    enu_expected_m = [ [ 0.0, 0.0, 100e6],
                       [ 0.0, 0.0, 100],
                       [ 0.0, 10e12 + p_rad, -m_rad ], # 0.0, distance + polar radius, meridonal radius

                       [ 0.0, p_rad, -m_rad ],
                       ] 
    
    # turn the test sites in to a list of topos objects
    topos_sites = [ Topos( latitude_degrees = _ts[0],
                           longitude_degrees = _ts[1],
                           elevation_m =_ts[2] ) for _ts in test_sites ]

    # turn the targets into a list of topos objects
    topos_targets = [ Topos( latitude_degrees = _tt[0],
                             longitude_degrees = _tt[1],
                             elevation_m = _tt[2] ) for _tt in test_targets ]

    now = datetime.now()
    now = now.replace(tzinfo=utc)
    
    # compute the east north up coordiates for each pair
    for i,site in enumerate( topos_sites ):
        _tt = topos_targets[i]
        topocentric = ( _tt - site).at( ts.utc( now) )
        enu_km = topocentric.east_north_up()
        enu_expect_km = np.asarray( enu_expected_m[i] ) / 1000.0
        tol = 100e-6 #  100 centimeter
        assert np.allclose( enu_km, enu_expect_km, rtol=0.0, atol=tol )
