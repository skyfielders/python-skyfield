import numpy as np
from skyfield import api
from skyfield.constants import DAY_S, tau
from skyfield.earthlib import earth_rotation_angle
from skyfield.framelib import true_equator_and_equinox_of_date
from skyfield.functions import A, from_spherical, length_of, mxv, rot_z
from skyfield.positionlib import Geocentric, ICRF, ITRF_to_GCRS2, _GIGAPARSEC_AU
from skyfield.starlib import Star
from .fixes import low_precision_ERA

from assay import assert_raises

def test_subtraction():
    p0 = ICRF((10,20,30), (40,50,60), center=0, target=499)
    p1 = ICRF((1,2,3), (4,5,6), center=0, target=399)
    p = p0 - p1
    assert p.center == 399
    assert p.target == 499
    assert isinstance(p, Geocentric)
    assert tuple(p.xyz.au) == (9, 18, 27)
    assert tuple(p.velocity.au_per_d) == (36, 45, 54)

    p1.center = 1
    with assert_raises(ValueError):
        p0 - p1

def test_separation_from_on_scalar():
    p0 = ICRF((1, 0, 0))
    p1 = ICRF((0, 1, 0))
    assert str(p0.separation_from(p1)) == '90deg 00\' 00.0"'

def test_separation_from_on_two_array_values():
    p0 = ICRF(([1,1], [0,0], [0,0]))
    p1 = ICRF(([0,-1], [1,0], [0,0]))
    sep = p0.separation_from(p1)
    d = sep.degrees
    assert len(d) == 2
    assert d[0] == 90.0
    assert d[1] == 180.0

def test_separation_from_on_an_array_and_a_scalar():
    p0 = ICRF(([1,0], [0,1], [0,0]))
    p1 = ICRF((0, 0, 1))
    sep = p0.separation_from(p1)
    d = sep.degrees
    assert len(d) == 2
    assert d[0] == 90.0
    assert d[1] == 90.0

    # And the other way around:

    sep = p1.separation_from(p0)
    d = sep.degrees
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

def test_frame_rotations_for_mean_of_date():
    ts = api.load.timescale()
    t = ts.utc(2020, 11, 21)
    p = ICRF((1.1,1.2,1.3), t=t)
    lat, lon, distance1 = p.frame_latlon(true_equator_and_equinox_of_date)

    # Verify that the frame_latlon() coordinates match those from the
    # more conventional radec() call.
    ra, dec, distance2 = p.radec(epoch='date')
    assert abs(lat.arcseconds() - dec.arcseconds()) < 1e-6
    assert abs(lon.arcseconds() - ra.arcseconds()) < 1e-6
    assert abs(distance1.au - distance2.au) < 1e-15

    # Now that we know the coordinates are good, we can use them to
    # rebuild a trusted x,y,z vector with which to test frame_xyz().
    x1, y1, z1 = from_spherical(distance1.au, lat.radians, lon.radians)
    x2, y2, z2 = p.frame_xyz(true_equator_and_equinox_of_date).au
    assert abs(x1 - x2) < 1e-15
    assert abs(y1 - y2) < 1e-15
    assert abs(z1 - z2) < 1e-15

def test_position_of_radec():
    epsilon = _GIGAPARSEC_AU * 1e-16

    p = api.position_of_radec(0, 0)
    assert length_of(p.xyz.au - [_GIGAPARSEC_AU, 0, 0]) < epsilon

    p = api.position_of_radec(6, 0)
    assert length_of(p.xyz.au - [0, _GIGAPARSEC_AU, 0]) < epsilon

    epsilon = 2e-16

    p = api.position_of_radec(12, 90, 2)
    assert length_of(p.xyz.au - [0, 0, 2]) < epsilon

    p = api.position_of_radec(12, 90, distance_au=2)
    assert length_of(p.xyz.au - [0, 0, 2]) < epsilon

    ts = api.load.timescale()
    epoch = ts.tt_jd(api.B1950)
    p = api.position_of_radec(0, 0, 1, epoch=epoch)
    assert length_of(p.xyz.au - [1, 0, 0]) > 1e-16
    ra, dec, distance = p.radec(epoch=epoch)
    assert abs(ra.hours) < 1e-12
    assert abs(dec.degrees) < 1e-12
    assert abs(distance.au - 1) < 3e-16

def test_position_from_radec():
    # Only a couple of minimal tests, since the routine is deprecated.
    p = api.position_from_radec(0, 0)
    assert length_of(p.xyz.au - [1, 0, 0]) < 1e-16

    p = api.position_from_radec(6, 0)
    assert length_of(p.xyz.au - [0, 1, 0]) < 1e-16

def test_ssb():
    ts = api.load.timescale()
    t = ts.utc(2025, 1, 28)
    p = api.SSB.at(t)
    z = [0,0,0]
    assert p.xyz.au.tolist() == z
    assert p.velocity.au_per_d.tolist() == z

    star = Star(ra_hours=12, dec_degrees=345)
    p.observe(star)

    t = ts.utc(2025, 1, [28,29])
    p = api.SSB.at(t)
    z2 = [[0,0], [0,0], [0,0]]
    assert p.xyz.au.tolist() == z2
    assert p.velocity.au_per_d.tolist() == z2

    p.observe(star)

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

    acceptable_error = 1e-11
    assert relative_error < acceptable_error

def test_light_time_method():
    p = ICRF([0.0, 1.0, 0.0])
    assert abs(p.light_time - 0.0057755183) < 1e-10

def test_hadec():
    # If the DE430 ephemeris excerpt is avaiable, this test can run
    # locally against the HA number from first line of
    # `moon_topo_4_6_2017_mkb_sf_v5_hadec.csv.txt` at:
    # https://github.com/skyfielders/python-skyfield/issues/510
    #planets = api.load('de430_1850-2150.bsp')
    #expected_ha = -0.660078756021

    # But in CI, we use DE421 for space and speed.
    planets = api.load('de421.bsp')
    expected_ha = -0.660078752

    ts = api.load.timescale()
    ts.polar_motion_table = [0.0], [0.009587], [0.384548]
    t = ts.utc(2017, 4, 6)
    topos = api.wgs84.latlon(-22.959748, -67.787260, elevation_m=5186.0)
    earth = planets['Earth']
    moon = planets['Moon']
    a = (earth + topos).at(t).observe(moon).apparent()
    ha, dec, distance = a.hadec()
    difference_mas = (ha.hours - expected_ha) * 15 * 3600 * 1e3
    assert abs(difference_mas) < 0.03

    # Drive-by test of position repr.
    assert repr(a) == (
        '<Apparent ICRS position and velocity at date t'
        ' center=WGS84 latitude -22.9597 N longitude -67.7873 E'
        ' elevation 5186.0 m target=301>'
    )

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
    assert np.allclose(tio_ra.degrees, era, rtol=0.0, atol=tol)
    assert np.allclose(tio_dec.degrees, 0.0, rtol=0.0, atol=tol)

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
    assert np.allclose(md_ra.degrees, era, rtol=0.0, atol=tol)
    assert np.allclose(md_dec.degrees, 90 - alt, rtol=0.0, atol=tol)

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

        assert np.allclose(ra_cirs.degrees, ra_sofa, rtol=0.0, atol=tol)
        assert np.allclose(dec_cirs.degrees, dec_sofa, rtol=0.0, atol=tol)

def test_phase_angle_and_fraction_illuminated():
    ts = api.load.timescale()
    t = ts.utc(2018, 9, range(9, 19), 5)
    t1 = t[-1]
    e = api.load('de421.bsp')
    earth, moon, sun = e['earth'], e['moon'], e['sun']

    p = earth.at(t1).observe(moon)
    a = p.phase_angle(sun).degrees.round(1)
    assert a == 76.2
    i = p.fraction_illuminated(sun).round(2)
    assert i == 0.62

    p = earth.at(t).observe(moon)
    a = p.phase_angle(sun).degrees.round(1)
    assert list(a) == [172.0, 172.6, 159.6, 146.6, 133.9,
                       121.7, 109.9, 98.4, 87.2, 76.2]
    i = (p.fraction_illuminated(sun) * 100).round(2)
    assert list(i) == [0.49, 0.41, 3.12, 8.25, 15.30, 23.73, 33.01,  # not 33.02?
                       42.71, 52.45, 61.92]

def test_astropy_conversion():
    try:
        import astropy
    except ImportError:
        # Drat: assay doesn't know about skipping a test.
        #raise SkipTest('AstroPy not installed')
        return
    else:
        astropy  # Use the library's name, to avoid a linter complaint.

    ts = api.load.timescale()
    r = np.array([1, 2, 3])
    t = ts.tt(2022, 1, 3)

    p = ICRF(r, t=t, center=0)
    a = p.to_skycoord()
    assert str(a) == '<SkyCoord (ICRS): (x, y, z) in AU\n    (1., 2., 3.)>'
    assert a.obstime is None

    p = ICRF(r, t=t, center=399)
    a = p.to_skycoord()
    assert str(a) == (
        '<SkyCoord (GCRS: obstime=2459582.5, obsgeoloc=(0., 0., 0.) m,'
        ' obsgeovel=(0., 0., 0.) m / s): (x, y, z) in AU\n    (1., 2., 3.)>'
    )
    assert a.obstime.fits == '2022-01-03T00:00:00.000'

    p = ICRF(r, t=t, center=3)
    with assert_raises(NotImplementedError):
        a = p.to_skycoord()

def test_old_position_attribute():
    ts = api.load.timescale()
    t = ts.tt(2022, 1, 3)
    r = A[6, 7, 8]
    p = ICRF(r, t=t, center=0)
    assert tuple(p.position.au) == (6, 7, 8)
