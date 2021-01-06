import os
from numpy import pi, seterr, linspace

from skyfield.api import load
from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN
from skyfield.data import mpc
from skyfield.keplerlib import _KeplerOrbit as KeplerOrbit, propagate, _CONVERT_GM
from skyfield.tests.test_elementslib import compare, ele_to_vec
from skyfield.units import Angle, Distance, Velocity

try:
    from io import BytesIO
except:
    from StringIO import StringIO as BytesIO

seterr(all='raise')

# Test against HORIZONS.

def test_against_horizons():
    # See the following files in the Skyfield repository:
    #
    # horizons/ceres-orbital-elements
    # horizons/ceres-position

    ts = load.timescale()
    t = ts.tdb_jd(2458886.500000000)

    a = 2.768873850275102E+00 # A
    e = 7.705857791518426E-02 # EC
    p_au = a * (1 - e*e)   # Wikipedia

    k = KeplerOrbit._from_mean_anomaly(
        semilatus_rectum_au=p_au,
        eccentricity=e,
        inclination_degrees=2.718528770987308E+01,
        longitude_of_ascending_node_degrees=2.336112629072238E+01,
        argument_of_perihelion_degrees=1.328964361683606E+02,
        mean_anomaly_degrees=1.382501360489816E+02,
        epoch=t,
        gm_km3_s2=GM_SUN,
        center=None,
        target=None,
    )
    r, v = k._at(t)[:2]
    sun_au = [
        -0.004105894975783999, 0.006739680703224941, 0.002956344702049446,
    ]
    horizons_au = [
        1.334875927366032E+00, -2.239607658161781E+00, -1.328895183461897E+00,
    ]
    epsilon = Distance(m=0.001).au
    assert abs(r + sun_au - horizons_au).max() < epsilon

def test_minor_planet():
    text = (b'00001    3.4   0.15 K205V 162.68631   73.73161   80.28698'
            b'   10.58862  0.0775571  0.21406009   2.7676569  0 MPO492748'
            b'  6751 115 1801-2019 0.60 M-v 30h Williams   0000      '
            b'(1) Ceres              20190915\n')

    ts = load.timescale()
    t = ts.utc(2020, 6, 17)
    eph = load('de421.bsp')
    df = mpc.load_mpcorb_dataframe(BytesIO(text))
    row = df.iloc[0]

    assert row.designation_packed == '00001'
    assert row.designation == '(1) Ceres'

    ceres = mpc.mpcorb_orbit(row, ts, GM_SUN)
    ra, dec, distance = eph['earth'].at(t).observe(eph['sun'] + ceres).radec()

    assert ceres.target == '(1) Ceres'
    assert abs(ra.hours - 23.1437) < 0.00005
    assert abs(dec.degrees - -17.323) < 0.0005

def test_comet():
    text = (b'    CJ95O010  1997 03 29.6333  0.916241  0.994928  130.6448'
            b'  283.3593   88.9908  20200224  -2.0  4.0  C/1995 O1 (Hale-Bopp)'
            b'                                    MPC106342\n')

    ts = load.timescale()
    t = ts.utc(2020, 5, 31)
    eph = load('de421.bsp')
    e = eph['earth'].at(t)

    for loader in mpc.load_comets_dataframe, mpc.load_comets_dataframe_slow:
        df = loader(BytesIO(text))
        row = df.iloc[0]
        k = mpc.comet_orbit(row, ts, GM_SUN)
        p = e.observe(eph['sun'] + k)
        ra, dec, distance = p.radec()

        # The file authorities/mpc-hale-bopp in the repository is the
        # source of these angles.  TODO: can we tighten this bound and
        # drive it to fractions of an arcsecond?

        ra_want = Angle(hours=(23, 59, 16.6))
        dec_want = Angle(degrees=(-84, 46, 58))
        assert abs(ra_want.arcseconds() - ra.arcseconds()) < 2.0
        assert abs(dec_want.arcseconds() - dec.arcseconds()) < 0.2
        assert abs(distance.au - 43.266) < 0.0005

        assert k.target == 'C/1995 O1 (Hale-Bopp)'

def test_comet_with_eccentricity_of_exactly_one():
    ts = load.timescale()
    t = ts.utc(2020, 8, 13)
    planets = load('de421.bsp')
    earth, sun = planets['earth'], planets['sun']

    data = (b'    CK15A020  2015 08  1.8353  5.341055  1.000000  208.8369  '
            b'258.5042  109.1696            10.5  4.0  C/2015 A2 (PANSTARRS)'
            b'                                    MPC 93587')

    with BytesIO(data) as f:
        df = mpc.load_comets_dataframe(f)

    df = df[df['designation'] == 'C/2015 A2 (PANSTARRS)']
    comet = mpc.comet_orbit(df.iloc[0], ts, GM_SUN)
    ra, dec, distance = earth.at(t).observe(sun + comet).radec()

    # These are exactly the RA and declination from the Minor Planet
    # Center for this comet!  (The RA seconds returned by Skyfield
    # actually say "46.45s", which would round up to 46.5, but what's a
    # tenth of an arcsecond among friends?)
    assert str(ra).startswith('18h 46m 46.4')
    assert str(dec).startswith("-72deg 05' 33.")

# Test various round-trips through the kepler orbit object.

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

def check_orbit(p, e, i, Om, w, v,
                p_eps=None, e_eps=None, i_eps=None, Om_eps=None, w_eps=None, v_eps=None):
    pos0, vel0 = ele_to_vec(p, e, i, Om, w, v, mu)

    pos1, vel1 = propagate(pos0, vel0, 0, times, mu)

    orbit = KeplerOrbit(Distance(km=pos1), Velocity(km_per_s=vel1), dummy_time, mu_au3_d2=mu*_CONVERT_GM)
    ele = orbit.elements_at_epoch

    if p_eps: compare(p, ele.semi_latus_rectum.km, p_eps)
    if e_eps: compare(e, ele.eccentricity, e_eps)
    if i_eps: compare(i, ele.inclination.radians, i_eps, mod=True)
    if Om_eps: compare(Om, ele.longitude_of_ascending_node.radians, Om_eps, mod=True)
    if w_eps: compare(w, ele.argument_of_periapsis.radians, w_eps, mod=True)
    if v_eps: compare(v, ele.true_anomaly.radians, v_eps, mod=True)


times = linspace(-1e11, 1e11, 1001) # -3170 years to +3170 years, including 0
mu = 403503.2355022598
dummy_time = load.timescale().utc(2018)

def test_circular():
    check_orbit(p=300000, e=0, i=.5, Om=1, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)


def test_circular_equatorial():
    check_orbit(p=300000, e=0, i=0, Om=0, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15)


def test_circular_retrograde_equatorial():
    check_orbit(p=300000, e=0, i=pi, Om=0, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15)


def test_circular_polar():
    check_orbit(p=300000, e=0, i=pi/2, Om=1, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)


def test_circular_non_zero_arg_of_periapsis():
    check_orbit(p=300000, e=0, i=.5, Om=1, w=.5, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)


def test_elliptical():
    check_orbit(p=300000, e=.3, i=1, Om=0, w=4, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)


def test_elliptical_equatorial():
    check_orbit(p=300000, e=.3, i=0, Om=0, w=1, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)


def test_elliptical_retrograde_equatorial():
    check_orbit(p=300000, e=.3, i=pi, Om=0, w=4, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)


def test_elliptical_polar():
    check_orbit(p=300000, e=.2, i=pi/2, Om=1, w=2, v=3,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-8)


def test_parabolic():
    check_orbit(p=300000, e=1, i=1, Om=0, w=4, v=3,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-13, Om_eps=1e-13, w_eps=1e-13)


def test_parabolic_equatorial():
    check_orbit(p=300000, e=1, i=0, Om=0, w=1, v=2,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-13)


def test_parabolic_retrograde_equatorial():
    check_orbit(p=300000, e=1, i=pi, Om=0, w=1, v=2,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-15, Om_eps=1e-13, w_eps=1e-13)


def test_parabolic_polar():
    check_orbit(p=300000, e=1, i=pi/2, Om=1, w=2, v=3,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-14, Om_eps=1e-13, w_eps=1e-13)


def test_hyperbolic():
    check_orbit(p=300000, e=1.3, i=1, Om=0, w=4, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-10, Om_eps=1e-10, w_eps=1e-6)


def test_hyperbolic_equatorial():
    check_orbit(p=300000, e=1.3, i=0, Om=0, w=1, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-6)


def test_hyperbolic_retrograde_equatorial():
    check_orbit(p=300000, e=1.3, i=pi, Om=0, w=1, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-15, Om_eps=1e-9, w_eps=1e-6)


def test_hyperbolic_polar():
    check_orbit(p=300000, e=1.3, i=pi/2, Om=1, w=2, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-10, Om_eps=1e-10, w_eps=1e-6)


def test_equatorial_non_zero_longitude_of_ascending_node():
    check_orbit(p=300000, e=.3, i=0, Om=0, w=4, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, w_eps=1e-7)
