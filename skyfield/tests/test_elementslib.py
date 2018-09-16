from skyfield.api import load, Time, load_file
from skyfield.data.spice import inertial_frames
from skyfield.units import Distance, Angle, Velocity
from skyfield.constants import DAY_S
from skyfield.elementslib import (
    OsculatingElements,
    normpi,
    osculating_elements_of,
)
from numpy import (array, sin, cos, pi, sqrt, ndarray, float64, repeat,
                   seterr, inf, linspace, arccos)
import os

def ts():
    yield load.timescale()

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

ECLIPTIC = inertial_frames['ECLIPJ2000']

ephem = load_file(_data_path('de430-2015-03-02.bsp'))
jup_ephem = load_file(_data_path('jup310-2015-03-02.bsp'))
io = jup_ephem['io']
moon = ephem['moon']
earth = ephem['earth']
sun = ephem['sun']

seterr(all='raise')


def compare(value, expected_value, epsilon, mod=False):
    """Compares value to expected value, and works if one or both are arrays.

    Also allows epsilon to be an array.

    If mod==True, then compare(0, tau, 0) is True.
    """
    if mod:
        diff = normpi(value - expected_value)
    else:
        diff = value - expected_value

    if hasattr(value, '__len__') or hasattr(expected_value, '__len__'):
        if hasattr(epsilon, '__len__'):
            assert (abs(diff) <= epsilon).all()
        else:
            assert max(abs(diff)) <= epsilon
    else:
        assert abs(diff) <= epsilon


def check_types(elements, length):
    """Checks that all of the attributes in the OsculatingElements object are
    present, have the correct type, and have the expected size.
    """
    for item in ['inclination', 'longitude_of_ascending_node',
                 'argument_of_periapsis', 'true_anomaly',
                 'argument_of_latitude', 'eccentric_anomaly', 'mean_anomaly',
                 'mean_longitude', 'true_longitude', 'longitude_of_periapsis',
                 'mean_motion_per_day']:
        element = getattr(elements, item)
        assert isinstance(element, Angle)
        assert element.radians.size == length

    for item in ['eccentricity', 'period_in_days']:
        element = getattr(elements, item)
        assert isinstance(element, (ndarray, float64))
        assert element.size == length

    for item in ['_n_vec', '_e_vec', '_h_vec']:
        element = getattr(elements, item)
        assert isinstance(element, (ndarray, float64))
        assert element[0].size == length

    for item in ['semi_latus_rectum', 'apoapsis_distance',
                 'periapsis_distance', 'semi_major_axis', 'semi_minor_axis']:
        element = getattr(elements, item)
        assert isinstance(element, Distance)
        assert element.au.size == length

        assert isinstance(elements.periapsis_time, Time)
        assert elements.periapsis_time.tt.size == length


def horizons_dict(elem, units='km_d', ):
    """
    Outputs dictionary with keys that match labels used by Horizons.

    The data from this method is exactly the same as the data in the
    elements object and differs from it only in units.

    Returns an dictionary whose keys match the labels from Horizons:
        EC = eccentricity
        QR = periapsis distance
        IN = inclination
        OM = longitude of ascending node
        W  = argument of periapsis
        Tp = time of closest periapsis
        N  = mean motion
        MA = mean anomaly
        TA = true anomaly
        A  = semi-major axis
        AD = apoapsis distance
        PR = period

    Just like in Horizons, all angle values are in degrees, but the
    distance and time units can be specified with the `units` keyword,
    and must be one of the following:
        'km_s' for kilometers and seconds
        'km_d' for kilometers and days
        'au_d' for au and days
    """
    data = {}
    data['EC'] = elem.eccentricity

    if units=='km_s' or units=='km_d':
        data['QR'] = elem.periapsis_distance.km
    elif units=='au_d':
        data['QR'] = elem.periapsis_distance.au

    data['IN'] = elem.inclination.degrees
    data['OM'] = elem.longitude_of_ascending_node.degrees
    data['W'] = elem.argument_of_periapsis.degrees
    data['Tp'] = elem.periapsis_time.tdb

    if units=='km_d' or units =='au_d':
        data['N'] = elem.mean_motion_per_day.degrees
    elif units=='km_s':
        data['N'] = elem.mean_motion_per_day.degrees/DAY_S

    data['MA'] = elem.mean_anomaly.degrees
    data['TA'] = elem.true_anomaly.degrees

    if units=='km_s' or units=='km_d':
        data['A'] = elem.semi_major_axis.km
        data['AD'] = elem.apoapsis_distance.km
    elif units=='au_d':
        data['A'] = elem.semi_major_axis.au
        data['AD'] = elem.apoapsis_distance.au

    if units=='km_d' or units =='au_d':
        data['PR'] = elem.period_in_days
    elif units=='km_s':
        data['PR'] = elem.period_in_days*DAY_S

    return data


def horizons_array(elem, units='km_d', ):
    """
    Outputs numpy array containing data in the same order as horizons.

    The data from this method is exactly the same as the data in the
    elements object and differs from it only in units.

    Just like in Horizons, all angle values are in degrees, but the
    distance and time units can be specified with the `units` keyword,
    and must be one of the following:
        'km_s' for kilometers and seconds
        'km_d' for kilometers and days
        'au_d' for au and days

    The shape of the array is ``(12,)`` if the time used to construct the
    position is a float, and ``(12, n)`` if the time is an array of length n.
    """
    dict_ = horizons_dict(elem, units=units)
    array_ = array([dict_['EC'],
                 dict_['QR'],
                 dict_['IN'],
                 dict_['OM'],
                 dict_['W'],
                 dict_['Tp'],
                 dict_['N'],
                 dict_['MA'],
                 dict_['TA'],
                 dict_['A'],
                 dict_['AD'],
                 dict_['PR']])
    return array_


def test_repr(ts):
    elements = osculating_elements_of((moon-earth).at(ts.utc(2015, 3, 2, 12)))
    assert repr(elements) == '<Elements 1 sets>'


def test_single_time(ts):
    """Tests creation of an OsculatingElements object with a single set of elements
    """
    elements = osculating_elements_of((moon-earth).at(ts.utc(2015, 3, 2, 12)))
    check_types(elements, 1)


def test_multiple_times(ts):
    """Tests creation of an OsculatingElements object with multiple sets of elements
    """
    time = ts.utc(2015, 3, 2, [12, 13, 14, 15])
    elements = osculating_elements_of((moon-earth).at(time))
    check_types(elements, len(time))


def test_equatorial_km_d(ts):
    """Tests against data from Horizons in km and days, with equatorial reference plane
    """
    geocentric_pos = (moon - earth).at(ts.tdb(2015, 3, 2, 2))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos))
    horizons_data = array([5.569337304355707E-02,
                           3.628019705296879E+05,
                           1.832507608006020E+01,
                           3.570946187642240E+02,
                           3.318968921151215E+02,
                           2457072.329517839011,
                           1.320459197700415E+01,
                           1.486020417866582E+02,
                           1.517384963232830E+02,
                           3.841993269696947E+05,
                           4.055966834097015E+05,
                           2.726324301628867E+01])
    epsilon = array([1e-10, 1e-4, 1e-7, 1e-6, 1e-6, 1e-8, 1e-9, 1e-6, 1e-7, 1e-4, 1e-4, 1e-8])
    compare(calc_data, horizons_data, epsilon)


def test_equatorial_km_s(ts):
    """Tests against data from Horizons in km and seconds, with equatorial reference plane
    """
    geocentric_pos = (moon - earth).at(ts.tdb(2015, 3, 2, 2))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos), units='km_s')
    horizons_data = array([5.569337304355707E-02,
                           3.628019705296879E+05,
                           1.832507608006020E+01,
                           3.570946187642240E+02,
                           3.318968921151215E+02,
                           2457072.329517839011,
                           1.528309256597703E-04,
                           1.486020417866582E+02,
                           1.517384963232830E+02,
                           3.841993269696947E+05,
                           4.055966834097015E+05,
                           2.355544196607342E+06])
    epsilon = array([1e-10, 1e-4, 1e-7, 1e-6, 1e-6, 1e-8, 1e-14, 1e-6, 1e-7, 1e-4, 1e-4, 1e-3])
    compare(calc_data, horizons_data, epsilon)


def test_equatorial_au_d(ts):
    """Tests against data from Horizons in au and days, with equatorial reference plane
    """
    geocentric_pos = (moon - earth).at(ts.tdb(2015, 3, 2, 2))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos), units='au_d')
    horizons_data = array([5.569337304355707E-02,
                           2.425181380136368E-03,
                           1.832507608006020E+01,
                           3.570946187642240E+02,
                           3.318968921151215E+02,
                           2457072.329517839011,
                           1.320459197700415E+01,
                           1.486020417866582E+02,
                           1.517384963232830E+02,
                           2.568213873445825E-03,
                           2.711246366755283E-03,
                           2.726324301628867E+01])
    epsilon = array([1e-10, 1e-12, 1e-7, 1e-6, 1e-6, 1e-8, 1e-9, 1e-6, 1e-7, 1e-13, 1e-12, 1e-8])
    compare(calc_data, horizons_data, epsilon)


def test_ecliptic_km_d(ts):
    """Tests against data from Horizons in km and days, with ecliptic reference plane
    """
    geocentric_pos = (moon - earth).at(ts.tdb(2015, 3, 2, 2))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos, ECLIPTIC), units='km_d')
    horizons_data = array([5.569337304355707E-02,
                          3.628019705296879E+05,
                          5.216521657765558E+00,
                          1.900948892867905E+02,
                          1.390846818575352E+02,
                          2457072.329517839011,
                          1.320459197700415E+01,
                          1.486020417866582E+02,
                          1.517384963232830E+02,
                          3.841993269696947E+05,
                          4.055966834097015E+05,
                          2.726324301628867E+01])
    epsilon = array([1e-10, 1e-4, 1e-8, 1e-6, 1e-6, 1e-8, 1e-9, 1e-6, 1e-7, 1e-4, 1e-4, 1e-8])
    compare(calc_data, horizons_data, epsilon)


def test_ecliptic_km_s(ts):
    """Tests against data from Horizons in km and seconds, with ecliptic reference plane
    """
    geocentric_pos = (moon - earth).at(ts.tdb(2015, 3, 2, 2))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos, ECLIPTIC), units='km_s')
    horizons_data = array([5.569337304355707E-02,
                           3.628019705296879E+05,
                           5.216521657765558E+00,
                           1.900948892867905E+02,
                           1.390846818575352E+02,
                           2457072.329517839011,
                           1.528309256597703E-04,
                           1.486020417866582E+02,
                           1.517384963232830E+02,
                           3.841993269696947E+05,
                           4.055966834097015E+05,
                           2.355544196607342E+06])
    epsilon = array([1e-10, 1e-4, 1e-8, 1e-6, 1e-6, 1e-8, 1e-14, 1e-6, 1e-7, 1e-4, 1e-4, 1e-3])
    compare(calc_data, horizons_data, epsilon)


def test_ecliptic_au_d(ts):
    """Tests against data from Horizons in au and days, with ecliptic reference plane
    """
    geocentric_pos = (moon - earth).at(ts.tdb(2015, 3, 2, 2))
    elem = osculating_elements_of(geocentric_pos, ECLIPTIC)
    calc_data = horizons_array(elem, units='au_d')
    horizons_data = array([5.569337304355707E-02,
                           2.425181380136368E-03,
                           5.216521657765558E+00,
                           1.900948892867905E+02,
                           1.390846818575352E+02,
                           2457072.329517839011,
                           1.320459197700415E+01,
                           1.486020417866582E+02,
                           1.517384963232830E+02,
                           2.568213873445825E-03,
                           2.711246366755283E-03,
                           2.726324301628867E+01])
    epsilon = array([1e-10, 1e-12, 1e-8, 1e-6, 1e-6, 1e-8, 1e-9, 1e-6, 1e-7, 1e-13, 1e-12, 1e-8])
    compare(calc_data, horizons_data, epsilon)


def test_extreme_ellipse(ts):
    """Tests against data from Horizons for an orbit with eccentricity just less than 1
    """
    geocentric_pos = (io - sun).at(ts.tdb(2015, 3, 2, 17, 26))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos), units='km_s')
    horizons_data = array([9.993434925710607E-01,
                           1.163126430217223E+08,
                           2.164979337400508E+01,
                           7.503612948248414E+00,
                           3.571949928607168E+02,
                           2456680.438054572791,
                           2.798965463884538E-10,
                           9.764838165348996E-03,
                           1.351769989470609E+02,
                           1.771688146920545E+11,
                           3.542213167410872E+11,
                           1.286189503390209E+12])
    epsilon = array([1e-10, 1e-1, 1e-8, 1e-8, 1e-7, 1e-8, 1e-16, 1e-8, 1e-7, 1e5, 1e5, 1e6])
    compare(calc_data, horizons_data, epsilon)


def test_slightly_hyperbolic(ts):
    """Tests against data from Horizons for an orbit with eccentricity just over 1
    """
    geocentric_pos = (io - sun).at(ts.tdb(2015, 3, 2, 17, 27))
    calc_data = horizons_array(osculating_elements_of(geocentric_pos), units='km_s')
    horizons_data = array([1.000249165282725E+00,
                           1.176022222580809E+08,
                           2.167088597088522E+01,
                           7.441033056578390E+00,
                           3.575776020085432E+02,
                           2456680.266815741546,
                           6.437043762782399E-11,
                           2.246667771669457E-03,
                           1.348525808471548E+02,
                           -4.719847844449893E+11])
    epsilon = array([1e-11, 1e-2, 1e-11, 1e-11, 1e-9, 1e-8, 1e-17, 1e-10, 1e-9, 1e4])
    compare(calc_data[:-2], horizons_data, epsilon)
    assert (calc_data[-2:] == array([inf, inf])).all()


def test_periapsis_time(ts):
    """This tests the moment when the eccentricity of Io orbiting the sun
    transitions from just below 1 to just above 1. Periapsis time is
    calculated using M/n, and both M and n go to 0 as e goes to 1.
    Periapsis time of Io around the sun should change linearly over a very
    small time interval interval (a fraction of one second) during this
    transition, but instead periapsis time deviated from linear within .005s
    of e=1.

    One source of the deviation was the fact that M was being remapped to
    -pi to pi when e>1. Using M without remapping reduced the deviation when
    e>1. The second thing that reduced the deviation was widening the
    tolerance for the use of the parabolic periapsis time equation.

    This test makes sure these fixes don't regress by checking that the
    periapsis times don't deviate from linear by more than 1e-5 days within
    .005s of e=1.
    """
    t = ts.tdb(jd=linspace(2457084.226893796, 2457084.226893912, 500))

    elem = osculating_elements_of((io - sun).at(t))
    Tp = elem.periapsis_time.tdb
    line = -240.61044176706827*t.tdb + 593656801.6052344
    compare(Tp, line, 1e-5)


def ele_to_vec(p, e, i, Om, w, v, mu):
    """Calculates state vectors from orbital elements. Also checks for invalid
    sets of elements.

    These equations are from this document:

    https://web.archive.org/web/*/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc

    """
    # Checks that longitude of ascending node is 0 if inclination is 0
    if isinstance(i, ndarray) or isinstance(Om, ndarray):
        if ((i==0)*(Om!=0)).any():
            raise ValueError('If inclination is 0, longitude of ascending node must be 0')
    else:
        if i==0 and Om!=0:
            raise ValueError('If inclination is 0, longitude of ascending node must be 0')

    # Checks that argument of periapsis is 0  if eccentricity is 0
    if isinstance(e, ndarray) or isinstance(w, ndarray):
        if ((e==0)*(w!=0)).any():
            raise ValueError('If eccentricity is 0, argument of periapsis must be 0')
    else:
        if e==0 and w!=0:
            raise ValueError('If eccentricity is 0, argument of periapsis must be 0')

    # Checks that true anomaly is less than arccos(-1/e) for hyperbolic orbits
    if isinstance(e, ndarray) and isinstance(v, ndarray):
        inds = (e>1)
        if (v[inds]>arccos(-1/e[inds])).any():
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')
    elif isinstance(e, ndarray) and not isinstance(v, ndarray):
        inds = (e>1)
        if (v>arccos(-1/e[inds])).any():
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')
    elif isinstance(v, ndarray) and not isinstance(e, ndarray):
        if e>1 and (v>arccos(-1/e)).any():
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')
    else:
        if e>1 and v>arccos(-1/e):
            raise ValueError('If eccentricity is >1, abs(true anomaly) cannot be more than arccos(-1/e)')

    # Checks that inclination is between 0 and pi
    if isinstance(i, ndarray):
        assert ((i>=0) * (i < pi)).all()
    else:
        assert i>=0 and i<pi
    r = p/(1 + e*cos(v))
    h = sqrt(p*mu)
    u = v+w

    X = r*(cos(Om)*cos(u) - sin(Om)*sin(u)*cos(i))
    Y = r*(sin(Om)*cos(u) + cos(Om)*sin(u)*cos(i))
    Z = r*(sin(i)*sin(u))

    X_dot = X*h*e/(r*p)*sin(v) - h/r*(cos(Om)*sin(u) + sin(Om)*cos(u)*cos(i))
    Y_dot = Y*h*e/(r*p)*sin(v) - h/r*(sin(Om)*sin(u) - cos(Om)*cos(u)*cos(i))
    Z_dot = Z*h*e/(r*p)*sin(v) + h/r*sin(i)*cos(u)

    # z and z_dot are independent of Om, so if Om is an array and the other
    # elements are scalars, z and z_dot need to be repeated
    if Z.size!=X.size:
        Z = repeat(Z, X.size)
        Z_dot = repeat(Z_dot, X.size)

    return array([X, Y, Z]), array([X_dot, Y_dot, Z_dot])


def check_orbit(p, e, i, Om, w, v, ts, mod=False):
    """Checks that the given set of elements are calculated properly by
    elementslib.py

    Converts the given elements to state vectors using ele_to_vec, then uses
    those state vectors to create an OsculatingElements object, and then
    checks that the data in the OsculatingElements object matches the input
    elements.
    """
    length = 1
    for item in p, e, i, Om, w, v:
        if isinstance(item, (list, ndarray)):
            length = len(item)

    mu = 403503.2355022598
    pos_vec, vel_vec = ele_to_vec(p, e, i, Om, w, v, mu)
    time_tt = ts.utc(2018).tt
    time = ts.tt(jd=repeat(time_tt, pos_vec[0].size))
    elements = OsculatingElements(Distance(km=pos_vec),
                                  Velocity(km_per_s=vel_vec),
                                  time, mu)
    check_types(elements, length)
    compare(time, elements.time, 1e-9)
    compare(elements.semi_latus_rectum.km, p, 1e-9)
    compare(elements.eccentricity, e, 1e-14)
    if mod:
        compare(elements.inclination.radians, i, 1e-14, mod=True)
        compare(elements.longitude_of_ascending_node.radians, Om, 1e-14, mod=True)
        compare(elements.argument_of_periapsis.radians, w, 1e-14, mod=True)
        compare(elements.true_anomaly.radians, v, 1e-14, mod=True)
    else:
        compare(elements.inclination.radians, i, 1e-14)
        compare(elements.longitude_of_ascending_node.radians, Om, 1e-14)
        compare(elements.argument_of_periapsis.radians, w, 1e-14)
        compare(elements.true_anomaly.radians, v, 1e-14)


# The remaining tests check that certain edge case orbits are being handled
# correctly. Each element is checked in the middle and at the edges of all
# quadrants.

angles1 = array([0, .25, .5, .75, 1, 1.25, 1.5, 1.75])*pi
angles2 = array([-.75, -.5, -.25, 0, .25, .5, .75])*pi

def test_circular(ts):
    e = 0
    w = 0
    check_orbit(300000, e, .5, angles1, w, 1, ts)
    check_orbit(300000, e, .5, 1, w, angles1, ts)

    for angle in angles1:
        check_orbit(300000, e, .5, angle, w, 1, ts)
        check_orbit(300000, e, .5, 1, w, angle, ts)

def test_circular_equatorial(ts):
    e = 0
    i = 0
    w = 0
    Om = 0
    check_orbit(300000, e, i, Om, w, angles1, ts)

    for angle in angles1:
        check_orbit(300000, e, i, Om, w, angle, ts)

def test_circular_polar(ts):
    e = 0
    w = 0
    i = pi/2
    check_orbit(300000, e, i, angles1, w, 3, ts)
    check_orbit(300000, e, i, 1, w, angles1, ts)

    for angle in angles1:
        check_orbit(300000, e, i, angle, w, 3, ts)
        check_orbit(300000, e, i, 1, w, angle, ts)

def test_elliptical(ts):
    check_orbit(300000, .3, angles1[:4], 0, 4, 5, ts, mod=True)
    check_orbit(300000, .3, .1, angles1, 4, 5, ts)
    check_orbit(300000, .3, .1, 2, angles1, 5, ts, mod=True)
    check_orbit(300000, .3, .1, 2, 4, angles1, ts)

    for angle in angles1[:4]:
        check_orbit(300000, .3, angle, 0, 4, 5, ts, mod=True)
    for angle in angles1:
        check_orbit(300000, .3, .1, angle, 4, 5, ts)
        check_orbit(300000, .3, .1, 2, angle, 5, ts, mod=True)
        check_orbit(300000, .3, .1, 2, 4, angle, ts)

def test_elliptical_equatorial(ts):
    i = 0
    Om = 0
    check_orbit(300000, .3, i, Om, angles1, 5, ts)
    check_orbit(300000, .3, i, Om, 4, angles1, ts, mod=True)

    for angle in angles1:
        check_orbit(300000, .3, i, Om, angle, 5, ts)
        check_orbit(300000, .3, i, Om, 4, angle, ts, mod=True)

def test_elliptical_polar(ts):
    i = pi/2
    check_orbit(300000, .2, i, angles1, 2, 3, ts, mod=True)
    check_orbit(300000, .2, i, 1, angles1, 3, ts, mod=True)
    check_orbit(300000, .2, i, 1, 2, angles1, ts)

    for angle in angles1:
        check_orbit(300000, .2, i, angle, 2, 3, ts, mod=True)
        check_orbit(300000, .2, i, 1, angle, 3, ts, mod=True)
        check_orbit(300000, .2, i, 1, 2, angle, ts)

def test_parabolic(ts):
    e = 1
    check_orbit(300000, e, angles1[:4], 0, 4, 3, ts)
    check_orbit(300000, e, 2, angles1, 4, 3, ts)
    check_orbit(300000, e, 2, 3, angles1, 3, ts)
    check_orbit(300000, e, 2, 3, 4, angles2, ts)

    for angle in angles1[:4]:
        check_orbit(300000, e, angle, 0, 4, 3, ts)
    for angle in angles1:
        check_orbit(300000, e, 2, angle, 4, 3, ts)
        check_orbit(300000, e, 2, 3, angle, 3, ts)
    for angle in angles2:
        check_orbit(300000, e, 2, 3, 4, angle, ts)

def test_parabolic_equatorial(ts):
    e = 1
    i = 0
    Om = 0
    check_orbit(300000, e, i, Om, angles1, 2, ts, mod=True)
    check_orbit(300000, e, i, Om, 4, angles2, ts)

    for angle in angles1:
        check_orbit(300000, e, i, Om, angle, 2, ts, mod=True)
    for angle in angles2:
        check_orbit(300000, e, i, Om, 4, angle, ts)

def test_parabolic_polar(ts):
    e = 1
    i = pi/2
    check_orbit(300000, e, i, angles1, 2, 3, ts)
    check_orbit(300000, e, i, 1, angles1, 3, ts)
    check_orbit(300000, e, i, 1, 2, angles2, ts)

    for angle in angles1:
        check_orbit(300000, e, i, angle, 2, 3, ts)
        check_orbit(300000, e, i, 1, angle, 3, ts)
    for angle in angles2:
        check_orbit(300000, e, i, 1, 2, angle, ts)

def test_hyperbolic(ts):
    e = 1.3
    check_orbit(300000, e, angles1[:4], 0, 4, .5, ts, mod=True)
    check_orbit(300000, e, 2, angles1, 4, .5, ts)
    check_orbit(300000, e, 2, 3, angles1, .5, ts)
    check_orbit(300000, e, 2, 3, 4, angles2, ts)

    for angle in angles1[:4]:
        check_orbit(300000, e, angle, 0, 4, .5, ts, mod=True)
    for angle in angles1:
        check_orbit(300000, e, 2, angle, 4, .5, ts)
        check_orbit(300000, e, 2, 3, angle, .5, ts)
    for angle in angles2:
        check_orbit(300000, e, 2, 3, 4, angle, ts)


def test_hyperbolic_equatorial(ts):
    e = 1.3
    i = 0
    Om = 0
    check_orbit(300000, e, i, Om, angles1, .5, ts)
    check_orbit(300000, e, i, Om, 4, angles2, ts)

    for angle in angles1:
        check_orbit(300000, e, i, Om, angle, .5, ts)
    for angle in angles2:
        check_orbit(300000, e, i, Om, 4, angle, ts)


def test_hyperbolic_polar(ts):
    e = 1.3
    i = pi/2
    check_orbit(300000, e, i, angles1, 2, .5, ts)
    check_orbit(300000, e, i, 1, angles1, .5, ts)
    check_orbit(300000, e, i, 1, 2, angles2, ts)

    for angle in angles1:
        check_orbit(300000, e, i, angle, 2, .5, ts)
        check_orbit(300000, e, i, 1, angle, .5, ts)
    for angle in angles2:
        check_orbit(300000, e, i, 1, 2, angle, ts)


def test_all_types_at_once(ts):

    check_orbit(p=array([300000]*12),
                e=array([ 0, 0,    0, .3, .3,   .2, 1, 1,    1, 1.3, 1.3,  1.3]),
                i=array([.5, 0, pi/2, .1,  0, pi/2, 2, 0, pi/2,   2,   0, pi/2]),
                Om=array([ 1, 0,    1,  2,  0,    1, 3, 0,    1,   3,   0,    1]),
                w=array([ 0, 0,    0,  4,  4,    2, 4, 4,    2,   4,   4,    2]),
                v=array([ 1, 5,    3,  5,  5,    3, 5, 5,    3,  .5,  .5,   .5]),
                ts=ts,
                mod=True)
