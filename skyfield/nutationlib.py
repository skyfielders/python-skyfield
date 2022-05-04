"""Routines that compute Earth nutation."""
from numpy import array, cos, dot, fmod, sin, outer, zeros
from .constants import ASEC2RAD, ASEC360, DEG2RAD, tau, T0
from .functions import load_bundled_npy

_TENTH_USEC_2_RAD = ASEC2RAD / 1e7
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

# These wrappers return nutation angles in radians as expected by the
# Time object.  We can't change the units returned by the underlying
# routines without breaking any applications that discovered them at
# some point in the past few years (though they are officially
# undocumented) and started calling them directly.

def iau2000a_radians(t, fundamental_argument_terms=5, lunisolar_terms=687,
                     planetary_terms=687):
    """Return the IAU 2000A angles delta-psi and delta-epsilon in radians."""
    d_psi, d_eps = iau2000a(t.tt, fundamental_argument_terms, lunisolar_terms,
                            planetary_terms)
    d_psi *= _TENTH_USEC_2_RAD
    d_eps *= _TENTH_USEC_2_RAD
    return d_psi, d_eps

def iau2000b_radians(t):
    """Return the IAU 2000B angles delta-psi and delta-epsilon in radians."""
    d_psi, d_eps = iau2000b(t.tt)
    d_psi *= _TENTH_USEC_2_RAD
    d_eps *= _TENTH_USEC_2_RAD
    return d_psi, d_eps

# Lower-level routines.

def build_nutation_matrix(mean_obliquity_radians,
                          true_obliquity_radians,
                          psi_radians):
    """Generate the nutation rotation matrix, given three nutation parameters.

    The input angles can be simple floats.  Or, they can be arrays of
    the same length, in which case the output matrix will have an extra
    dimension of that same length providing *n* rotation matrices.

    """
    cobm = cos(mean_obliquity_radians)
    sobm = sin(mean_obliquity_radians)
    cobt = cos(true_obliquity_radians)
    sobt = sin(true_obliquity_radians)
    cpsi = cos(psi_radians)
    spsi = sin(psi_radians)

    return array(((cpsi,
                  -spsi * cobm,
                  -spsi * sobm),
                  (spsi * cobt,
                   cpsi * cobm * cobt + sobm * sobt,
                   cpsi * sobm * cobt - cobm * sobt),
                  (spsi * sobt,
                   cpsi * cobm * sobt - sobm * cobt,
                   cpsi * sobm * sobt + cobm * cobt)))

def mean_obliquity(jd_tdb):
    """Return the mean obliquity of the ecliptic in arcseconds.

    The caller need only supply a single argument:

    `jd_tdb` - TDB time as a Julian date float, or NumPy array of floats

    The formulae used to compute the mean obliquity are based on
    equations 37 and 39 from:

    Capitaine et al. (2003), _Astronomy and Astrophysics_ 412, 567-586.

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

    This routine takes a single argument:

    `jd_tt` - Terrestrial Time: Julian date float, or NumPy array of floats

    The formulae used are from:

    Capitaine, N., Wallace, P.T., and McCarthy, D.D. (2003). _Astron. &
    Astrophys._ 406, p. 1135-1149. Table 3.

    _IERS Conventions (2010)_, Chapter 5, p. 60, Table 5.2e.  (Table
    5.2e presented in the printed publication is a truncated series. The
    full series, which is used here, is available on the IERS
    Conventions Center website in file tab5.2e.txt.)
    ftp://tai.bipm.org/iers/conv2010/chapter5/

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

    a = ke1.dot(fa)
    c_terms = se1_0 * sin(a)
    c_terms += se1_1 * cos(a)
    c_terms *= t

    a = ke0_t.dot(fa)
    c_terms += se0_t_0.dot(sin(a))
    c_terms += se0_t_1.dot(cos(a))

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

def iau2000a(jd_tt, fundamental_argument_terms=5, lunisolar_terms=687,
             planetary_terms=687):
    """Compute Earth nutation based on the IAU 2000A nutation model.

    ``jd_tt`` - Terrestrial Time: Julian date float, or NumPy array of floats

    Returns a tuple ``(delta_psi, delta_epsilon)`` measured in tenths of
    a micro-arcsecond.  Each value is either a float, or a NumPy array
    with the same dimensions as the input argument.

    Supply smaller integer values for ``fundamental_argument_terms``,
    ``lunisolar_terms``, and ``planetary_terms`` to trade off accuraccy
    for speed.

    """
    # Interval between fundamental epoch J2000.0 and given date.

    t = (jd_tt - T0) / 36525.0

    # Compute fundamental arguments from Simon et al. (1994), in radians.

    a = fundamental_arguments(t, fundamental_argument_terms)

    # ** Luni-solar nutation **
    # Summation of luni-solar nutation series.

    cutoff = lunisolar_terms
    arg = nals_t[:cutoff].dot(a).T

    sarg = sin(arg)
    carg = cos(arg)

    dpsi = dot(sarg, lunisolar_longitude_coefficients[:cutoff,0])
    dpsi += dot(sarg, lunisolar_longitude_coefficients[:cutoff,1]) * t
    dpsi += dot(carg, lunisolar_longitude_coefficients[:cutoff,2])

    deps = dot(carg, lunisolar_obliquity_coefficients[:cutoff,0])
    deps += dot(carg, lunisolar_obliquity_coefficients[:cutoff,1]) * t
    deps += dot(sarg, lunisolar_obliquity_coefficients[:cutoff,2])

    # Compute and add in planetary components.

    if not planetary_terms:
        return dpsi, deps

    if getattr(t, 'shape', ()) == ():
        a = t * anomaly_coefficient + anomaly_constant
    else:
        a = (outer(anomaly_coefficient, t).T + anomaly_constant).T
    a[-1] *= t

    cutoff = planetary_terms
    arg = napl_t[:cutoff].dot(a).T

    sarg = sin(arg)
    carg = cos(arg)

    dpsi += dot(sarg, nutation_coefficients_longitude[:cutoff,0])
    dpsi += dot(carg, nutation_coefficients_longitude[:cutoff,1])

    deps += dot(sarg, nutation_coefficients_obliquity[:cutoff,0])
    deps += dot(carg, nutation_coefficients_obliquity[:cutoff,1])

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
    dpsi, deps = iau2000a(jd_tt, 2, 77, 0)
    dpsi += -0.000135e7
    deps +=  0.000388e7
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

def fundamental_arguments(t, terms=5):
    """Compute the fundamental arguments (mean elements) of Sun and Moon.

    ``t`` - TDB time in Julian centuries since J2000.0, as float or NumPy array

    Outputs fundamental arguments, in radians:
          a[0] = l (mean anomaly of the Moon)
          a[1] = l' (mean anomaly of the Sun)
          a[2] = F (mean argument of the latitude of the Moon)
          a[3] = D (mean elongation of the Moon from the Sun)
          a[4] = Omega (mean longitude of the Moon's ascending node);
                 from Simon section 3.4(b.3),
                 precession = 5028.8200 arcsec/cy)

    Pass a smaller value for the number of polynomial ``terms`` if you
    want to trade accuracy for speed.

    """
    fa = iter((fa4, fa3, fa2, fa1)[-terms+1:])
    a = next(fa) * t
    for fa_i in fa:
        a += fa_i
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

# Deprecated functions that third-party code might still call; several
# of our tests also still call them, to help keep them working.

def compute_nutation(t):
    """Deprecated: this is now a method on the Time object."""
    return t.N

def earth_tilt(t):
    """Deprecated: these are now computed separately on the Time object."""
    d_psi, d_eps = t._nutation_angles_radians
    mean_ob = t._mean_obliquity_radians
    true_ob = mean_ob + d_eps
    c_terms = equation_of_the_equinoxes_complimentary_terms(t.tt)
    eq_eq = d_psi * cos(mean_ob) + c_terms
    return (mean_ob / DEG2RAD, true_ob / DEG2RAD, eq_eq / ASEC2RAD / 15.0,
            d_psi / ASEC2RAD, d_eps / ASEC2RAD)
