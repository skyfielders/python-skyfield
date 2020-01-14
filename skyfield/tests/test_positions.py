import numpy as np
from skyfield import api, constants
from skyfield.earthlib import earth_rotation_angle
from skyfield.functions import length_of
from skyfield.positionlib import ICRF
from skyfield.starlib import Star

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

def test_position_from_radec():
    p = api.position_from_radec(0, 0)
    assert length_of(p.position.au - [1, 0, 0]) < 1e-16

    p = api.position_from_radec(6, 0)
    assert length_of(p.position.au - [0, 1, 0]) < 1e-16

    p = api.position_from_radec(12, 90, 2)
    assert length_of(p.position.au - [0, 0, 2]) < 2e-16

    ts = api.load.timescale(builtin=True)
    epoch = ts.tt_jd(api.B1950)
    p = api.position_from_radec(0, 0, epoch=epoch)
    assert length_of(p.position.au - [1, 0, 0]) > 1e-16
    ra, dec, distance = p.radec(epoch=epoch)
    assert abs(ra.hours) < 1e-12
    assert abs(dec.degrees) < 1e-12
    assert abs(distance.au - 1) < 1e-16

def test_itrf_vector():
    ts = api.load.timescale(builtin=True)
    t = ts.utc(2019, 11, 2, 3, 53)
    top = api.Topos(latitude_degrees=45, longitude_degrees=0,
                    elevation_m=constants.AU_M - constants.ERAD)
    x, y, z = top.at(t).itrf_xyz().au
    assert abs(x - 0.7071) < 1e-4
    assert abs(y - 0.0) < 1e-14
    assert abs(z - 0.7071) < 1e-4

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
