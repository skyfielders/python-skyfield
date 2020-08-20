from __future__ import division

import sys
import math
from numpy import(abs, amax, amin, arange, arccos, arctan, array, cos, cosh,
                  cross, exp, log, ndarray, newaxis, ones_like, pi, power,
                  repeat, sin, sinh, sqrt, sum, tan, tanh, tile, zeros_like)

from skyfield.constants import AU_KM, DAY_S, DEG2RAD
from skyfield.functions import dots, length_of, mxv
from skyfield.descriptorlib import reify
from skyfield.elementslib import OsculatingElements, normpi
from skyfield.units import Distance, Velocity
from skyfield.vectorlib import VectorFunction

_CONVERT_GM = DAY_S * DAY_S / AU_KM / AU_KM / AU_KM

class _KeplerOrbit(VectorFunction):
    def __init__(self,
                 position,
                 velocity,
                 epoch,
                 mu_au3_d2,
                 center=None,
                 target=None,
        ):
        """ Calculates the position of an object using 2 body propagation

        Parameters
        ----------
        position : Distance
            Position vector at epoch with shape (3,)
        velocity : Velocity
            Velocity vector at epoch with shape (3,)
        epoch : Time
            Time corresponding to `position` and `velocity`
        mu_au_d : float
            Value of mu (G * M) in au^3/d^2
        center : int
            NAIF ID of the primary body, 399 for geocentric orbits, 10 for
            heliocentric orbits
        target : int
            NAIF ID of the secondary body
        """
        self.position_at_epoch = position
        self.velocity_at_epoch = velocity
        self.epoch = epoch
        self.mu_au3_d2 = mu_au3_d2
        self.center = center
        self.target = target

        self._rotation = None  # TODO: make argument?

    @classmethod
    def _from_periapsis(
            cls,
            semilatus_rectum_au,
            eccentricity,
            inclination_degrees,
            longitude_of_ascending_node_degrees,
            argument_of_perihelion_degrees,
            t_periapsis,
            gm_km3_s2,
            center=None,
            target=None,
        ):
        """Build a `KeplerOrbit` given its parameters and date of periapsis."""
        gm_au3_d2 = gm_km3_s2 * _CONVERT_GM
        pos, vel = ele_to_vec(
            semilatus_rectum_au,
            eccentricity,
            DEG2RAD * inclination_degrees,
            DEG2RAD * longitude_of_ascending_node_degrees,
            DEG2RAD * argument_of_perihelion_degrees,
            0.0,
            gm_au3_d2,
        )
        return cls(
            Distance(pos),
            Velocity(vel),
            t_periapsis,
            gm_au3_d2,
            center,
            target,
        )

    @classmethod
    def _from_true_anomaly(cls, p, e, i, Om, w, v,
                          epoch,
                          mu_km_s=None,
                          mu_au3_d2=None,
                          center=None,
                          target=None,
        ):
        """ Creates a `KeplerOrbit` object from elements using true anomaly

        Parameters
        ----------
        p : Distance
            Semi-Latus Rectum
        e : float
             Eccentricity
        i : Angle
            Inclination
        Om : Angle
            Longitude of Ascending Node
        w : Angle
            Argument of periapsis
        v : Angle
            True anomaly
        epoch : Time
            Time corresponding to `position` and `velocity`
        mu_km_s : float
            Value of mu (G * M) in km^3/s^2
        mu_au3_d2 : float
            Value of mu (G * M) in au^3/d^2
        center : int
            NAIF ID of the primary body, 399 for geocentric orbits, 10 for
            heliocentric orbits
        target : int
            NAIF ID of the secondary body
        """
        if (mu_km_s and mu_au3_d2) or (not mu_km_s and not mu_au3_d2):
            raise ValueError('Either mu_km_s or mu_au3_d2 should be used, but not both')

        if mu_au3_d2:
            mu_km_s = mu_au3_d2 * AU_KM**3 / DAY_S**2

        position, velocity = ele_to_vec(p.km,
                                        e,
                                        i.radians,
                                        Om.radians,
                                        w.radians,
                                        v.radians,
                                        mu_km_s,
        )
        return cls(Distance(km=position),
                   Velocity(km_per_s=velocity),
                   epoch,
                   mu_km_s,
                   center=center,
                   target=target,
        )


    @classmethod
    def _from_mean_anomaly(
            cls,
            semilatus_rectum_au,
            eccentricity,
            inclination_degrees,
            longitude_of_ascending_node_degrees,
            argument_of_perihelion_degrees,
            mean_anomaly_degrees,
            epoch,
            gm_km3_s2,
            center=None,
            target=None,
        ):
        """ Creates a `KeplerOrbit` object from elements using mean anomaly

        Parameters
        ----------
        p : Distance
            Semi-Latus Rectum
        e : float
             Eccentricity
        i : Angle
            Inclination
        Om : Angle
            Longitude of Ascending Node
        w : Angle
            Argument of periapsis
        M : Angle
            Mean anomaly
        epoch : Time
            Time corresponding to `position` and `velocity`
        mu_km_s : float
            Value of mu (G * M) in km^3/s^2
        mu_au3_d2 : float
            Value of mu (G * M) in au^3/d^2
        center : int
            NAIF ID of the primary body, 399 for geocentric orbits, 10 for
            heliocentric orbits
        target : int
            NAIF ID of the secondary body
        """
        M = DEG2RAD * mean_anomaly_degrees
        gm_au3_d2 = gm_km3_s2 * _CONVERT_GM
        if eccentricity < 1.0:
            E = eccentric_anomaly(eccentricity, M)
            v = true_anomaly_closed(eccentricity, E)
        elif eccentricity > 1.0:
            E = eccentric_anomaly(eccentricity, M)
            v = true_anomaly_hyperbolic(eccentricity, E)
        else:
            v = true_anomaly_parabolic(semilatus_rectum_au, gm_au3_d2, M)

        pos, vel = ele_to_vec(
            semilatus_rectum_au,
            eccentricity,
            DEG2RAD * inclination_degrees,
            DEG2RAD * longitude_of_ascending_node_degrees,
            DEG2RAD * argument_of_perihelion_degrees,
            v,
            gm_au3_d2,
        )
        return cls(
            Distance(pos),
            Velocity(vel),
            epoch,
            gm_au3_d2,
            center,
            target,
        )

    def _at(self, time):
        """Propagate the KeplerOrbit to the given Time object

        The Time object can contain one time, or an array of times
        """
        pos, vel = propagate(
            self.position_at_epoch.au,
            self.velocity_at_epoch.au_per_d,
            self.epoch.tt,
            time.tt,
            self.mu_au3_d2,
        )
        if self._rotation is not None:
            pos = mxv(self._rotation, pos)
            vel = mxv(self._rotation, vel)
        return pos, vel, None, None


    @reify
    def elements_at_epoch(self):
        return OsculatingElements(self.position_at_epoch,
                                  self.velocity_at_epoch,
                                  self.epoch,
                                  mu_km_s = self.mu_au3_d2 / _CONVERT_GM,
        )


    def __str__(self):
        ele = self.elements_at_epoch
        if self.target_name:
            return 'KeplerOrbit {0} {1} -> {2} {3}'.format(self.center,
                                                           self.center_name,
                                                           self.target,
                                                           self.target_name,
            )
        else:
            ele = self.elements_at_epoch
            string = 'KeplerOrbit {0} {1} -> q={2:.2}au e={3:.3f} i={4:.1f} Om={5:.1f} w={6:.1f}'
            return string.format(self.center,
                                 self.center_name,
                                 ele.periapsis_distance.au,
                                 ele.eccentricity,
                                 ele.inclination.degrees,
                                 ele.longitude_of_ascending_node.degrees,
                                 ele.argument_of_periapsis.degrees,
            )


    def __repr__(self):
        return '<{0}>'.format(str(self))


def eccentric_anomaly(e, M):
    """ Iterates to solve Kepler's equation to find eccentric anomaly

    Based on the algorithm in section 8.10.2 of the Explanatory Supplement
    to the Astronomical Almanac, 3rd ed.
    """
    M = normpi(M)
    E = M + e*sin(M)

    max_iters = 100
    iters = 0
    while iters < max_iters:
        dM = M - (E - e*sin(E))
        dE = dM/(1 - e*cos(E))
        E = E + dE
        if abs(dE) < 1e-14: return E
        iters += 1
    else:
        raise ValueError('Failed to converge')


def true_anomaly_hyperbolic(e, E):
    """Calculates true anomaly from eccentricity and eccentric anomaly.

    Valid for hyperbolic orbits. Equations from the relevant Wikipedia entries.

    """
    return 2.0 * arctan(sqrt((e + 1.0) / (e - 1.0)) * tanh(E/2))


def true_anomaly_closed(e, E):
    """Calculates true anomaly from eccentricity and eccentric anomaly.

    Valid for closed orbits. Equations from the relevant Wikipedia entries.

    """
    return 2.0 * arctan(sqrt((1.0 + e) / (1.0 - e)) * tan(E/2))


def true_anomaly_parabolic(p, gm, M):
    """Calculates true anomaly from semi-latus rectum, gm, and mean anomaly.

    Valid for parabolic orbits. Equations from
    https://en.wikipedia.org/wiki/Parabolic_trajectory.

    """
    delta_t = sqrt(2 * p**3 / gm) * M # from http://www.bogan.ca/orbits/kepler/orbteqtn.html
    periapsis_distance = p / 2
    A = 3 / 2 * sqrt(gm / (2 * periapsis_distance**3)) * delta_t
    B = (A + (A**2 + 1))**(1/3)
    return 2 * arctan(B - 1/B)


def ele_to_vec(p, e, i, Om, w, v, mu):
    """Calculates state vectors from orbital elements. Also checks for invalid
    sets of elements.

    Based on equations from this document:

    https://web.archive.org/web/*/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    """
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

    # Checks that inclination is in the range [0, pi]
    if isinstance(i, ndarray):
        if not ((i>=0) * (i <= pi)).all():
            raise ValueError('Inclination outside the range [0, pi] radians')
    else:
        if not 0 <= i <= pi:
            raise ValueError('Inclination outside the range [0, pi] radians')

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


dpmax = sys.float_info.max

def bracket(num, end1, end2):
    num[num<end1] = end1
    num[num>end2] = end2
    return num


def find_trunc():
    denom = 2
    factr = 2
    trunc = 1
    x = 1 / denom
    while 1+x > 1:
        denom = denom * (2+factr) * (1+factr)
        factr = factr + 2
        trunc = trunc + 1
        x = 1 / denom
    return trunc


trunc = find_trunc()
odd_factorials = array([math.factorial(i) for i in range(3, trunc*2, 2)])
even_factorials = array([math.factorial(i) for i in range(2, trunc*2, 2)])


def stumpff(x):
    """Calculates Stumpff functions

    Based on the function toolkit/src/spicelib/stmp03.f from the SPICE toolkit,
    which can be downloaded from naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
    """
    if (x < (-(log(2) + log(dpmax))**2)).any():
        raise ValueError('Argument below lower bound')

    z = sqrt(abs(x))

    c0 = zeros_like(x)
    c1 = zeros_like(x)
    c2 = zeros_like(x)
    c3 = zeros_like(x)

    low = x < -1
    c0[low] = cosh(z[low])
    c1[low] = sinh(z[low])/z[low]
    c2[low] = (1 - c0[low])/x[low]
    c3[low] = (1 - c1[low])/x[low]

    high = x > 1
    c0[high] = cos(z[high])
    c1[high] = sin(z[high])/z[high]
    c2[high] = (1 - c0[high])/x[high]
    c3[high] = (1 - c1[high])/x[high]

    mid = ~low * ~high
    n = sum(mid)
    exponents = tile(arange(0, trunc-1), [n, 1])
    odd_denominators = tile(odd_factorials, [n, 1])
    even_denominators = tile(even_factorials, [n, 1])
    numerators = repeat(x[mid][newaxis].T, trunc-1, axis=1)
    c3[mid] = (sum(power(numerators[:, ::2], exponents[:, ::2])/odd_denominators[:, ::2], axis=1)
        - sum(power(numerators[:, 1::2], exponents[:, 1::2])/odd_denominators[:, 1::2], axis=1))
    c2[mid] = (sum(power(numerators[:, ::2], exponents[:, ::2])/even_denominators[:, ::2], axis=1)
        - sum(power(numerators[:, 1::2], exponents[:, 1::2])/even_denominators[:, 1::2], axis=1))
    c1[mid] = 1 - x[mid]*c3[mid]
    c0[mid] = 1 - x[mid]*c2[mid]

    return c0, c1, c2, c3


def propagate(position, velocity, t0, t1, gm):
    """Propagates a position and velocity vector with an array of times.

    Based on the function toolkit/src/spicelib/prop2b.f from the SPICE toolkit,
    which can be downloaded from naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html

    Parameters
    ----------
    position : ndarray
       Position vector with shape (3,)
    velocity : ndarray
        Velocity vector with shape (3,)
    t0 : float
        Time corresponding to `position` and `velocity`
    t1 : float or ndarray
        Time or times to propagate to
    gm : float
        Gravitational parameter in units that match the other arguments
    """
    if gm <= 0:
        raise ValueError("'gm' should be positive")
    if length_of(velocity) == 0:
        raise ValueError('Velocity vector has zero magnitude')
    if length_of(position) == 0:
        raise ValueError('Position vector has zero magnitude')

    r0 = length_of(position)
    rv = dots(position, velocity)

    hvec = cross(position, velocity)
    h2 = dots(hvec, hvec)

    if h2 == 0:
        raise ValueError('Motion is not conical')

    eqvec = cross(velocity, hvec)/gm + -position/r0
    e = length_of(eqvec)
    q = h2 / (gm * (1+e))

    f = 1 - e
    b = sqrt(q/gm)

    br0 = b * r0
    b2rv = b**2 * rv
    bq = b * q
    qovr0 = q / r0


    maxc = max(abs(br0),
               abs(b2rv),
               abs(bq),
               abs(qovr0/(bq)))

    if f < 0:
        fixed = log(dpmax/2) - log(maxc)
        rootf = sqrt(-f)
        logf = log(-f)
        bound = min(fixed/rootf, (fixed + 1.5*logf)/rootf)
    else:
        logbound = (log(1.5) + log(dpmax) - log(maxc)) / 3
        bound = exp(logbound)

    def kepler(x):
        c0, c1, c2, c3 = stumpff(f*x*x)
        return x*(br0*c1 + x*(b2rv*c2 + x*(bq*c3)))

    dt = t1 - t0

    if not isinstance(dt, ndarray):
        dt = array([dt])
        return_1d_array = True
    else:
        return_1d_array = False

    x = bracket(dt/bq, -bound, bound)
    kfun = kepler(x)

    past = dt < 0
    future = dt > 0

    upper = zeros_like(dt, dtype='float64')
    lower = zeros_like(dt, dtype='float64')
    oldx = zeros_like(dt, dtype='float64')

    lower[past] = x[past]
    upper[future] = x[future]

    while (kfun[past] > dt[past]).any():
        upper[past] = lower[past]
        lower[past] *= 2
        oldx[past] = x[past]
        x[past] = bracket(lower[past], -bound, bound)
        if (x[past] == oldx[past]).any():
            raise ValueError('The input delta time (dt) has a value of {0}.'
                             'This is beyond the range of DT for which we '
                             'can reliably propagate states. The limits for '
                             'this GM and initial state are from {1}'
                             'to {2}.'.format(dt, kepler(-bound), kepler(bound)))
        kfun[past] = kepler(x[past])

    while (kfun[future] < dt[future]).any():
        lower[future] = upper[future]
        upper[future] *= 2
        oldx[future] = x[future]
        x[future] = bracket(upper[future], -bound, bound)
        if (x[future] == oldx[future]).any():
            raise ValueError('The input delta time (dt) has a value of {0}.'
                             'This is beyond the range of DT for which we '
                             'can reliably propagate states. The limits for '
                             'this GM and initial state are from {1} '
                             'to {2}.'.format(dt, kepler(-bound), kepler(bound)))
        kfun[future] = kepler(x[future])

    x = amin(array([upper, amax(array([lower, (lower+upper)/2]), axis=0)]), axis=0)

    lcount = zeros_like(dt)
    mostc = ones_like(dt)*1000
    not_done = (lower < x) * (x < upper)

    while not_done.any():
        kfun[not_done] = kepler(x[not_done])

        high = (kfun > dt) * not_done
        low = (kfun < dt) * not_done
        same = (~high * ~low) * not_done

        upper[high] = x[high]
        lower[low] = x[low]
        upper[same] = lower[same] = x[same]

        condition = not_done * (mostc > 64) * (upper != 0) * (lower != 0)
        mostc[condition] = 64
        lcount[condition] = 0

        # vectorized version of min(upper, max(lower, (upper + lower)/2))
        x[not_done] = amin(array([upper[not_done], amax(array([lower[not_done], (lower[not_done]+upper[not_done])/2]), axis=0)]), axis=0)
        lcount += 1
        not_done = (lower < x) * (x < upper) * (lcount < mostc)

    c0, c1, c2, c3 = stumpff(f*x*x)
    br = br0*c0 + x*(b2rv*c1 + x*(bq*c2))

    pc = 1 - qovr0 * x**2 * c2
    vc = dt - bq * x**3 * c3
    pcdot = -qovr0 / br * x * c1
    vcdot = 1 - bq / br * x**2 * c2

    if return_1d_array:
        position_prop = pc*position + vc*velocity
        velocity_prop = pcdot*position + vcdot*velocity
    else:
        position_prop = pc*tile(position[newaxis].T, dt.size) + vc*tile(velocity[newaxis].T, dt.size)
        velocity_prop = pcdot*tile(position[newaxis].T, dt.size) + vcdot*tile(velocity[newaxis].T, dt.size)

    return position_prop, velocity_prop
