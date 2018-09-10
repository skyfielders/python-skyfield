"""Routines that compute Earth nutation."""
from numpy import array, cos, fmod, sin, outer, tensordot, zeros
from .constants import ASEC2RAD, ASEC360, DEG2RAD, tau, T0
from .functions import load_bundled_npy

_arrays = load_bundled_npy('nutation.npz')

ke0_t = _arrays['ke0_t']
ke1 = _arrays['ke1']
lunisolar_longitude_coefficients = _arrays['lunisolar_longitude_coefficients']
lunisolar_obliquity_coefficients = _arrays['lunisolar_obliquity_coefficients']
nals_t = _arrays['nals_t']
napl_t = _arrays['napl_t']
nutation_coefficients_longitude = _arrays['nutation_coefficients_longitude']
nutation_coefficients_obliquity = _arrays['nutation_coefficients_obliquity']
se0_t_0 = _arrays['se0_t_0']
se0_t_1 = _arrays['se0_t_1']

def compute_nutation(t):
    """Generate the nutation rotations for Time `t`.

    If the Julian date is scalar, a simple ``(3, 3)`` matrix is
    returned; if the date is an array of length ``n``, then an array of
    matrices is returned with dimensions ``(3, 3, n)``.

    """
    oblm, oblt, eqeq, psi, eps = t._earth_tilt

    cobm = cos(oblm * DEG2RAD)
    sobm = sin(oblm * DEG2RAD)
    cobt = cos(oblt * DEG2RAD)
    sobt = sin(oblt * DEG2RAD)
    cpsi = cos(psi * ASEC2RAD)
    spsi = sin(psi * ASEC2RAD)

    return array(((cpsi,
                  -spsi * cobm,
                  -spsi * sobm),
                  (spsi * cobt,
                   cpsi * cobm * cobt + sobm * sobt,
                   cpsi * sobm * cobt - cobm * sobt),
                  (spsi * sobt,
                   cpsi * cobm * sobt - sobm * cobt,
                   cpsi * sobm * sobt + cobm * cobt)))

def earth_tilt(t):
    """Return a tuple of information about the earth's axis and position.

    `t` - A Time object.

    The returned tuple contains five items:

    ``mean_ob`` - Mean obliquity of the ecliptic in degrees.
    ``true_ob`` - True obliquity of the ecliptic in degrees.
    ``eq_eq`` - Equation of the equinoxes in seconds of time.
    ``d_psi`` - Nutation in longitude in arcseconds.
    ``d_eps`` - Nutation in obliquity in arcseconds.

    """
    dp, de = t._nutation_angles
    c_terms = equation_of_the_equinoxes_complimentary_terms(t.tt) / ASEC2RAD

    d_psi = dp * 1e-7 + t.psi_correction
    d_eps = de * 1e-7 + t.eps_correction

    mean_ob = mean_obliquity(t.tdb)
    true_ob = mean_ob + d_eps

    mean_ob /= 3600.0
    true_ob /= 3600.0

    eq_eq = d_psi * cos(mean_ob * DEG2RAD) + c_terms
    eq_eq /= 15.0

    return mean_ob, true_ob, eq_eq, d_psi, d_eps

#

def mean_obliquity(jd_tdb):
    """Return the mean obliquity of the ecliptic in arcseconds.

    `jd_tt` - TDB time as a Julian date float, or NumPy array of floats

    """
    # Compute time in Julian centuries from epoch J2000.0.

    t = (jd_tdb - T0) / 36525.0

    # Compute the mean obliquity in arcseconds.  Use expression from the
    # reference's eq. (39) with obliquity at J2000.0 taken from eq. (37)
    # or Table 8.

    epsilon = (((( -  0.0000000434   * t
                   -  0.000000576  ) * t
                   +  0.00200340   ) * t
                   -  0.0001831    ) * t
                   - 46.836769     ) * t + 84381.406

    return epsilon

def equation_of_the_equinoxes_complimentary_terms(jd_tt):
    """Compute the complementary terms of the equation of the equinoxes.

    `jd_tt` - Terrestrial Time: Julian date float, or NumPy array of floats

    """
    # Interval between fundamental epoch J2000.0 and current date.

    t = (jd_tt - T0) / 36525.0

    # Build array for intermediate results.

    shape = getattr(jd_tt, 'shape', ())
    fa = zeros((14,) if shape == () else (14, shape[0]))

    # Mean Anomaly of the Moon.

    fa[0] = ((485868.249036 +
              (715923.2178 +
              (    31.8792 +
              (     0.051635 +
              (    -0.00024470)
              * t) * t) * t) * t) * ASEC2RAD
              + (1325.0*t % 1.0) * tau)

    # Mean Anomaly of the Sun.

    fa[1] = ((1287104.793048 +
              (1292581.0481 +
              (     -0.5532 +
              (     +0.000136 +
              (     -0.00001149)
              * t) * t) * t) * t) * ASEC2RAD
              + (99.0*t % 1.0) * tau)

    # Mean Longitude of the Moon minus Mean Longitude of the Ascending
    # Node of the Moon.

    fa[2] = (( 335779.526232 +
              ( 295262.8478 +
              (    -12.7512 +
              (     -0.001037 +
              (      0.00000417)
              * t) * t) * t) * t) * ASEC2RAD
              + (1342.0*t % 1.0) * tau)

    # Mean Elongation of the Moon from the Sun.

    fa[3] = ((1072260.703692 +
              (1105601.2090 +
              (     -6.3706 +
              (      0.006593 +
              (     -0.00003169)
              * t) * t) * t) * t) * ASEC2RAD
              + (1236.0*t % 1.0) * tau)

    # Mean Longitude of the Ascending Node of the Moon.

    fa[4] = (( 450160.398036 +
              (-482890.5431 +
              (      7.4722 +
              (      0.007702 +
              (     -0.00005939)
              * t) * t) * t) * t) * ASEC2RAD
              + (-5.0*t % 1.0) * tau)

    fa[ 5] = (4.402608842 + 2608.7903141574 * t)
    fa[ 6] = (3.176146697 + 1021.3285546211 * t)
    fa[ 7] = (1.753470314 +  628.3075849991 * t)
    fa[ 8] = (6.203480913 +  334.0612426700 * t)
    fa[ 9] = (0.599546497 +   52.9690962641 * t)
    fa[10] = (0.874016757 +   21.3299104960 * t)
    fa[11] = (5.481293872 +    7.4781598567 * t)
    fa[12] = (5.311886287 +    3.8133035638 * t)
    fa[13] = (0.024381750 +    0.00000538691 * t) * t

    fa %= tau

    # Evaluate the complementary terms.

    a = ke0_t.dot(fa)
    s0 = se0_t_0.dot(sin(a)) + se0_t_1.dot(cos(a))

    a = ke1.dot(fa)
    s1 = se1_0 * sin(a) + se1_1 * cos(a)

    c_terms = s0 + s1 * t
    c_terms *= ASEC2RAD
    return c_terms

anomaly_constant, anomaly_coefficient = array([

    # Mean anomaly of the Moon.
    (2.35555598, 8328.6914269554),

    # Mean anomaly of the Sun.
    (6.24006013, 628.301955),

    # Mean argument of the latitude of the Moon.
    (1.627905234, 8433.466158131),

    # Mean elongation of the Moon from the Sun.
    (5.198466741, 7771.3771468121),

    # Mean longitude of the ascending node of the Moon.
    (2.18243920, - 33.757045),

    # Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
    (4.402608842, 2608.7903141574),
    (3.176146697, 1021.3285546211),
    (1.753470314,  628.3075849991),
    (6.203480913,  334.0612426700),
    (0.599546497,   52.9690962641),
    (0.874016757,   21.3299104960),
    (5.481293871,    7.4781598567),
    (5.321159000,    3.8127774000),

    # General accumulated precession in longitude (gets multiplied by t).
    (0.02438175, 0.00000538691),
    ]).T

def iau2000a(jd_tt):
    """Compute Earth nutation based on the IAU 2000A nutation model.

    `jd_tt` - Terrestrial Time: Julian date float, or NumPy array of floats

    Returns a tuple ``(delta_psi, delta_epsilon)`` measured in tenths of
    a micro-arcsecond.  Each value is either a float, or a NumPy array
    with the same dimensions as the input argument.

    """
    # Interval between fundamental epoch J2000.0 and given date.

    t = (jd_tt - T0) / 36525.0

    # Compute fundamental arguments from Simon et al. (1994), in radians.

    a = fundamental_arguments(t)

    # ** Luni-solar nutation **
    # Summation of luni-solar nutation series (in reverse order).

    arg = nals_t.dot(a)
    fmod(arg, tau, out=arg)

    sarg = sin(arg)
    carg = cos(arg)

    stsc = array((sarg, t * sarg, carg)).T
    ctcs = array((carg, t * carg, sarg)).T

    dpsi = tensordot(stsc, lunisolar_longitude_coefficients)
    deps = tensordot(ctcs, lunisolar_obliquity_coefficients)

    # Compute and add in planetary components.

    if getattr(t, 'shape', ()) == ():
        a = t * anomaly_coefficient + anomaly_constant
    else:
        a = (outer(anomaly_coefficient, t).T + anomaly_constant).T
    a[-1] *= t

    fmod(a, tau, out=a)
    arg = napl_t.dot(a)
    fmod(arg, tau, out=arg)
    sc = array((sin(arg), cos(arg))).T

    dpsi += tensordot(sc, nutation_coefficients_longitude)
    deps += tensordot(sc, nutation_coefficients_obliquity)

    return dpsi, deps

def iau2000b(jd_tt):
    """Compute Earth nutation based on the faster IAU 2000B nutation model.

    `jd_tt` - Terrestrial Time: Julian date float, or NumPy array of floats

    Returns a tuple ``(delta_psi, delta_epsilon)`` measured in tenths of
    a micro-arcsecond.  Each is either a float, or a NumPy array with
    the same dimensions as the input argument.  The result will not take
    as long to compute as the full IAU 2000A series, but should still
    agree with ``iau2000a()`` to within a milliarcsecond between the
    years 1995 and 2020.

    """
    dpplan = -0.000135 * 1e7
    deplan =  0.000388 * 1e7

    t = (jd_tt - T0) / 36525.0

    # TODO: can these be replaced with fa0 and f1?

    el  = fmod (485868.249036 +
          t * 1717915923.2178, ASEC360) * ASEC2RAD;

    elp = fmod (1287104.79305 +
             t * 129596581.0481, ASEC360) * ASEC2RAD;

    f   = fmod (335779.526232 +
             t * 1739527262.8478, ASEC360) * ASEC2RAD;

    d   = fmod (1072260.70369 +
             t * 1602961601.2090, ASEC360) * ASEC2RAD;

    om  = fmod (450160.398036 -
             t * 6962890.5431, ASEC360) * ASEC2RAD;

    a = array((el, elp, f, d, om))

    arg = nals_t[:77].dot(a)
    fmod(arg, tau, out=arg)

    sarg = sin(arg)
    carg = cos(arg)

    stsc = array((sarg, t * sarg, carg)).T
    ctcs = array((carg, t * carg, sarg)).T

    dp = tensordot(stsc, lunisolar_longitude_coefficients[:77,])
    de = tensordot(ctcs, lunisolar_obliquity_coefficients[:77,])

    dpsi = dpplan + dp
    deps = deplan + de

    return dpsi, deps

fa0, fa1, fa2, fa3, fa4 = array((

    # Mean Anomaly of the Moon.
    (485868.249036, 1717915923.2178, 31.8792, 0.051635, - .00024470),

    # Mean Anomaly of the Sun.
    (1287104.79305,  129596581.0481, - 0.5532, 0.000136, - 0.00001149),

    # Mean Longitude of the Moon minus Mean Longitude of the Ascending
    # Node of the Moon.
    (335779.526232, 1739527262.8478, - 12.7512, -  0.001037, 0.00000417),

    # Mean Elongation of the Moon from the Sun.
    (1072260.70369, 1602961601.2090, - 6.3706, 0.006593, - 0.00003169),

    # Mean Longitude of the Ascending Node of the Moon.
    (450160.398036, - 6962890.5431, 7.4722, 0.007702, - 0.00005939),

    )).T[:,:,None]

def fundamental_arguments(t):

    """Compute the fundamental arguments (mean elements) of Sun and Moon.

    `t` - TDB time in Julian centuries since J2000.0, as float or NumPy array

    Outputs fundamental arguments, in radians:
          a[0] = l (mean anomaly of the Moon)
          a[1] = l' (mean anomaly of the Sun)
          a[2] = F (mean argument of the latitude of the Moon)
          a[3] = D (mean elongation of the Moon from the Sun)
          a[4] = Omega (mean longitude of the Moon's ascending node);
                 from Simon section 3.4(b.3),
                 precession = 5028.8200 arcsec/cy)

    """
    a = fa4 * t
    a += fa3
    a *= t
    a += fa2
    a *= t
    a += fa1
    a *= t
    a += fa0
    fmod(a, ASEC360, out=a)
    a *= ASEC2RAD
    if getattr(t, 'shape', ()):
        return a
    return a[:,0]

# Sine and cosine coefficients for t^1.

se1_0 = -0.87e-6
se1_1 = +0.00e-6
