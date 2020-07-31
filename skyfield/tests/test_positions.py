import numpy as np
from skyfield import api
from skyfield.constants import DAY_S, tau
from skyfield.earthlib import earth_rotation_angle
from skyfield.functions import length_of, mxv, rot_z
from skyfield.positionlib import ICRF, ITRF_to_GCRS2, _GIGAPARSEC_AU
from skyfield.starlib import Star
from .fixes import IS_32_BIT, low_precision_ERA

def test_subtraction():
    p0 = ICRF((10,20,30), (40,50,60))
    p1 = ICRF((1,2,3), (4,5,6))
    p = p0 - p1
    assert tuple(p.position.au) == (9, 18, 27)
    assert tuple(p.velocity.au_per_d) == (36, 45, 54)

def test_separation_from_on_scalar():
    p0 = ICRF((1, 0, 0))
    p1 = ICRF((0, 1, 0))
    assert str(p0.separation_from(p1)) == '90deg 00\' 00.0"'

def test_separation_from_on_two_array_values():
    p0 = ICRF(([1,1], [0,0], [0,0]))
    p1 = ICRF(([0,-1], [1,0], [0,0]))
    sep = p0.separation_from(p1)
    d = sep._degrees
    assert len(d) == 2
    assert d[0] == 90.0
    assert d[1] == 180.0

def test_separation_from_on_an_array_and_a_scalar():
    p0 = ICRF(([1,0], [0,1], [0,0]))
    p1 = ICRF((0, 0, 1))
    sep = p0.separation_from(p1)
    d = sep._degrees
    assert len(d) == 2
    assert d[0] == 90.0
    assert d[1] == 90.0

    # And the other way around:

    sep = p1.separation_from(p0)
    d = sep._degrees
    assert len(d) == 2
    assert d[0] == 90.0
    assert d[1] == 90.0

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

def test_position_of_radec():
    epsilon = _GIGAPARSEC_AU * 1e-16

    p = api.position_of_radec(0, 0)
    assert length_of(p.position.au - [_GIGAPARSEC_AU, 0, 0]) < epsilon

    p = api.position_of_radec(6, 0)
    assert length_of(p.position.au - [0, _GIGAPARSEC_AU, 0]) < epsilon

    epsilon = 2e-16

    p = api.position_of_radec(12, 90, 2)
    assert length_of(p.position.au - [0, 0, 2]) < epsilon

    p = api.position_of_radec(12, 90, distance_au=2)
    assert length_of(p.position.au - [0, 0, 2]) < epsilon

    ts = api.load.timescale()
    epoch = ts.tt_jd(api.B1950)
    p = api.position_of_radec(0, 0, 1, epoch=epoch)
    assert length_of(p.position.au - [1, 0, 0]) > 1e-16
    ra, dec, distance = p.radec(epoch=epoch)
    assert abs(ra.hours) < 1e-12
    assert abs(dec.degrees) < 1e-12
    assert abs(distance.au - 1) < 1e-16

def test_position_from_radec():
    # Only a couple of minimal tests, since the routine is deprecated.
    p = api.position_from_radec(0, 0)
    assert length_of(p.position.au - [1, 0, 0]) < 1e-16

    p = api.position_from_radec(6, 0)
    assert length_of(p.position.au - [0, 1, 0]) < 1e-16

def test_velocity_in_ITRF_to_GCRS2():
    # TODO: Get test working with these vectors too, showing it works
    # with a non-zero velocity vector, but in that case the test will
    # have to be fancier in how it corrects.
    # r = np.array([(1, 0, 0), (1, 1 / DAY_S, 0)]).T
    # v = np.array([(0, 1, 0), (0, 1, 0)]).T

    ts = api.load.timescale()
    t = ts.utc(2020, 7, 17, 8, 51, [0, 1])
    r = np.array([(1, 0, 0), (1, 0, 0)]).T
    v = np.array([(0, 0, 0), (0, 0, 0)]).T

    r, v = ITRF_to_GCRS2(t, r, v, True)

    # Rotate back to equinox-of-date before applying correction.
    r = mxv(t.M, r)
    v = mxv(t.M, v)

    r0, r1 = r.T
    v0 = v[:,0]

    # Apply a correction: the instantaneous velocity does not in fact
    # carry the position in a straight line, but in an arc around the
    # origin; so use trigonometry to move the destination point to where
    # linear motion would have carried it.
    angvel = (t.gast[1] - t.gast[0]) / 24.0 * tau
    r1 = mxv(rot_z(np.arctan(angvel) - angvel), r1)
    r1 *= np.sqrt(1 + angvel*angvel)

    actual_motion = r1 - r0
    predicted_motion = v0 / DAY_S

    relative_error = (length_of(actual_motion - predicted_motion)
                      / length_of(actual_motion))

    acceptable_error = 4e-12 if IS_32_BIT else 2e-12
    assert relative_error < acceptable_error

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
        with low_precision_ERA():
            ra_cirs, dec_cirs, _ = earth.at(st).observe(ss).apparent().cirs_radec(st)

        assert np.allclose(ra_cirs._degrees, ra_sofa, rtol=0.0, atol=tol)
        assert np.allclose(dec_cirs._degrees, dec_sofa, rtol=0.0, atol=tol)
