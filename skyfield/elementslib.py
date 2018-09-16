"""Computation of osculating elements that matches NASA HORIZONS."""

from .functions import dots, length_of, angle_between
from .constants import DAY_S, tau
from .data.gravitational_parameters import GM_dict
from .units import Distance, Angle, Velocity
from .descriptorlib import reify
from numpy import (array, arctan2, sin, arctan, tan, inf, repeat, float64,
                   sinh, sqrt, arccos, arctanh, zeros_like, ones_like, divide,
                   where, pi, cross)

def osculating_elements_of(position, reference_frame=None):
    """Produce the osculating orbital elements for a position.

    The ``position`` should be an :class:`~skyfield.positionlib.ICRF`
    instance like that returned by the ``at()`` method of any Solar
    System body, specifying a position, a velocity, and a time.  An
    instance of :class:`~skyfield.elementslib.OsculatingElements` is
    returned.

    """
    mu = GM_dict.get(position.center, 0) + GM_dict.get(position.target, 0)
    
    if reference_frame is not None:
        position_vec = Distance(reference_frame.dot(position.position.au))
        velocity_vec = Velocity(reference_frame.dot(position.velocity.au_per_d))
    else:
        position_vec = position.position
        velocity_vec = position.velocity
        
    return OsculatingElements(position_vec,
                              velocity_vec,
                              position.t,
                              mu)


class OsculatingElements(object):
    """
    Contains one or more sets of osculating orbital elements. Different
    elements are accessed as attributes of the ``OsculatingElements`` object.

    An ``OsculatingElements`` object can be initialized with the following
    parameters:

    position : Distance object
        Position vector with shape (3,) or (3, n)
    velocity : Velocity object
        Velocity vector with shape (3,) or (3, n)
    time: Time object
        The times of the position and velocity vectors
    mu_km_s: float
        Gravitational parameter (G*M) in units of km^3/s^2

    """
    def __init__(self, position, velocity, time, mu_km_s):
        self._pos_vec = position.km
        self._vel_vec = velocity.km_per_s
        self.time = time
        self._mu = mu_km_s

        self._h_vec = cross(self._pos_vec, self._vel_vec, 0, 0).T
        self._e_vec = eccentricity_vector(self._pos_vec, self._vel_vec, self._mu)
        self._n_vec = node_vector(self._h_vec)

    @reify
    def apoapsis_distance(self):
        Q = apoapsis_distance(self.semi_latus_rectum.km, self.eccentricity)
        return Distance(km=Q)

    @reify
    def argument_of_latitude(self):
        w = self.argument_of_periapsis.radians
        v = self.true_anomaly.radians
        u = (w+v)%tau
        return Angle(radians=u)

    @reify
    def argument_of_periapsis(self):
        w = argument_of_periapsis(self._n_vec, self._e_vec, self._pos_vec, self._vel_vec)
        return Angle(radians=w)

    @reify
    def eccentric_anomaly(self):
        E = eccentric_anomaly(self.true_anomaly.radians,
                              self.eccentricity,
                              self.semi_latus_rectum.km)
        return Angle(radians=E)

    @reify
    def eccentricity(self):
        return length_of(self._e_vec)

    @reify
    def inclination(self):
        i = inclination(self._h_vec)
        return Angle(radians=i)

    @reify
    def longitude_of_ascending_node(self):
        Om = longitude_of_ascending_node(self.inclination.radians, self._h_vec)
        return Angle(radians=Om)

    @reify
    def longitude_of_periapsis(self):
        Om = self.longitude_of_ascending_node.radians
        w = self.argument_of_periapsis.radians
        lp = (Om + w)%tau
        return Angle(radians=lp)

    @reify
    def mean_anomaly(self):
        M = mean_anomaly(self.eccentric_anomaly.radians, self.eccentricity)
        return Angle(radians=M)

    @reify
    def mean_longitude(self):
        L = (self.longitude_of_ascending_node.radians
              + self.argument_of_periapsis.radians
              + self.mean_anomaly.radians)%tau
        return Angle(radians=L)

    @reify
    def mean_motion_per_day(self):
        n = mean_motion(self.semi_major_axis.km,
                        self._mu)
        return Angle(radians=n*DAY_S)

    @reify
    def periapsis_distance(self):
        q = periapsis_distance(self.semi_latus_rectum.km, self.eccentricity)
        return Distance(km=q)

    @reify
    def periapsis_time(self):
        M = mean_anomaly(self.eccentric_anomaly.radians, self.eccentricity, shift=False)
        tp = time_since_periapsis(M,
                                  self.mean_motion_per_day.radians/DAY_S,
                                  self.true_anomaly.radians,
                                  self.semi_latus_rectum.km,
                                  self._mu)
        ts = self.time.ts
        times = self.time.tdb - tp/DAY_S
        return ts.tdb(jd=times)

    @reify
    def period_in_days(self):
        P = period(self.semi_major_axis.km, self._mu)
        return P/DAY_S

    @reify
    def semi_latus_rectum(self):
        p = semi_latus_rectum(self._h_vec, self._mu)
        return Distance(km=p)

    @reify
    def semi_major_axis(self):
        a = semi_major_axis(self.semi_latus_rectum.km, self.eccentricity)
        return Distance(km=a)

    @reify
    def semi_minor_axis(self):
        b = semi_minor_axis(self.semi_latus_rectum.km, self.eccentricity)
        return Distance(km=b)

    @reify
    def true_anomaly(self):
        v = true_anomaly(self._e_vec, self._pos_vec, self._vel_vec, self._n_vec)
        return Angle(radians=v)

    @reify
    def true_longitude(self):
        Om = self.longitude_of_ascending_node.radians
        w = self.argument_of_periapsis.radians
        v = self.true_anomaly.radians
        l = (Om + w + v)%tau
        return Angle(radians=l)

    def __repr__(self):
        return '<Elements {0} sets>'.format(self.time.tt.size)

# a = semi-major axis
# b = semi-minor axis
# e = eccentricity
# E = eccentric anomaly
# h = specific angular momentum
# i = inclination
# l = true longitude
# L = mean longitude
# M = mean anomaly
# n = mean motion
# Om = longitude of ascending node
# p = semi-latus rectum
# P = period
# q = periapsis distance
# Q = apoapsis distance
# t = time
# u = argument of latitude
# v = true anomaly
# w = argument of periapsis
# lp = longitude of periapsis

# Sources:
# Mostly this book:
# Bate, Mueller, & White. Fundamentals of Astrodynamics (1971), Section 2.4, pgs. 61-64

# Also these sites:
# https://web.archive.org/web/*/http://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
# https://web.archive.org/web/*/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
# https://en.wikipedia.org/wiki/Orbital_elements
# http://www.bogan.ca/orbits/kepler/orbteqtn.html


def normpi(num):
    return (num + pi)%tau - pi


def apoapsis_distance(p, e):
    if p.ndim == 0:
        return p*(1+e)/(1-e**2) if e < 1 else float64(inf)
    else:
        return divide(p*(1+e), 1-e**2, out=repeat(inf, p.shape), where=e<(1-1e-15))


def argument_of_periapsis(n_vec, e_vec, pos_vec, vel_vec):
    if n_vec.ndim == 1:
        if length_of(e_vec) < 1e-15: # circular
            return 0

        elif length_of(n_vec) < 1e-15: # equatorial and not circular
            angle = arctan2(e_vec[1], e_vec[0])%tau
            return angle if cross(pos_vec, vel_vec, 0, 0).T[2] >= 0 else -angle%tau

        else: # not circular and not equatorial
            angle = angle_between(n_vec, e_vec)
            return angle if e_vec[2] > 0 else -angle%tau
    else:
        w = zeros_like(pos_vec[0]) # defaults to 0 for circular orbits

        equatorial = length_of(n_vec) < 1e-15
        circular = length_of(e_vec) < 1e-15

        inds = ~circular*equatorial
        angle = arctan2(e_vec[1][inds], e_vec[0][inds])%tau
        condition = cross(pos_vec[:, inds], vel_vec[:, inds], 0, 0).T[2] >= 0
        w[inds] = where(condition, angle, -angle%tau)

        inds = ~circular*~equatorial
        angle = angle_between(n_vec[:, inds], e_vec[:, inds])
        condition = e_vec[2][inds] > 0
        w[inds] = where(condition, angle, -angle%tau)
        return w


def eccentric_anomaly(v, e, p):
    if e.ndim == 0:
        if e < 1:
            return 2*arctan(sqrt((1-e)/(1+e)) * tan(v/2))
        elif e > 1:
            return normpi(2*arctanh(tan(v/2)/sqrt((e+1)/(e-1))))
        else:
            return 0
    else:
        E = zeros_like(e) # defaults to 0 for parabolic

        inds = (e<1) # elliptical
        E[inds] = 2*arctan(sqrt((1-e[inds])/(1+e[inds])) * tan(v[inds]/2))

        inds = (e>1) # hyperbolic
        E[inds] = normpi(2*arctanh(tan(v[inds]/2)/sqrt((e[inds]+1)/(e[inds]-1))))

        return E


def eccentricity(h, a, mu):
    condition = (h**2/(a*mu) <= 1)
    if h.ndim == 0:
        return sqrt(1 - h**2/(a*mu)) if condition else float64(0)
    else:
        return sqrt(1 - h**2/(a*mu), out=zeros_like(h), where=condition)


def eccentricity_vector(pos_vec, vel_vec, mu):
    r = length_of(pos_vec)
    v = length_of(vel_vec)
    return ((v**2 - mu/r)*pos_vec - dots(pos_vec, vel_vec)*vel_vec)/mu


def inclination(h_vec):
    k_vec = array([zeros_like(h_vec[0]), zeros_like(h_vec[0]), ones_like(h_vec[0])])
    return angle_between(h_vec, k_vec)


def longitude_of_ascending_node(i, h_vec):
    if i.ndim == 0:
        return arctan2(h_vec[0], -h_vec[1])%tau if i != 0 else float64(0)
    else:
        return arctan2(h_vec[0], -h_vec[1], out=zeros_like(i), where=i!=0)%tau


def mean_anomaly(E, e, shift=True):
    if e.ndim == 0:
        if e < 1:
            return (E - e*sin(E))%tau
        elif e > 1:
            M = e*sinh(E) - E
            return normpi(M) if shift else M
        else:
            return float64(0)
    else:
        M = zeros_like(e) # defaults to 0 for parabolic

        inds = (e<1) # elliptical
        M[inds] = (E[inds] - e[inds]*sin(E[inds]))%tau

        inds = (e>1) # hyperbolic
        if shift:
            M[inds] = normpi(e[inds]*sinh(E[inds]) - E[inds])
        else:
            M[inds] = e[inds]*sinh(E[inds]) - E[inds]

        return M


def mean_motion(a, mu):
    return sqrt(mu/abs(a)**3)


def node_vector(h_vec):
    n_vec = array([-h_vec[1], h_vec[0], zeros_like(h_vec[0])]) # h_vec cross [0, 0, 1]
    n = length_of(n_vec)

    if h_vec.ndim == 1:
        return n_vec/n if n!=0 else n_vec
    else:
        return divide(n_vec, n, out=n_vec, where=n!=0)


def periapsis_distance(p, e):
    if p.ndim == 0:
        return p*(1-e) / (1-e**2) if e!=1 else p/2
    else:
        return divide(p*(1-e), (1-e**2), out=p/2, where=e!=1)


def period(a, mu):
    if a.ndim == 0:
        return tau*sqrt(a**3/mu) if a>0 else float64(inf)
    else:
        return tau*sqrt(a**3/mu, out=repeat(inf, a.shape), where=a>0)


def semi_latus_rectum(h_vec, mu):
    return length_of(h_vec)**2/mu


def semi_major_axis(p, e):
    if p.ndim == 0:
        return p/(1 - e**2) if e!=1 else float64(inf)
    else:
        return divide(p, 1 - e**2, out=repeat(inf, p.shape), where=e!=1)


def semi_minor_axis(p, e):
    if e.ndim == 0:
        if e < 1:
            return p/sqrt(1 - e**2)
        elif e > 1:
            return p*sqrt(e**2 - 1) / (1 - e**2)
        else:
            return float64(0)
    else:
        b = zeros_like(e) # defaults to 0 for parabolic

        inds = (e<1) # elliptical
        b[inds] = p[inds]/sqrt(1 - e[inds]**2)

        inds = (e>1) # hyperbolic
        b[inds] = p[inds]*sqrt(e[inds]**2 - 1) / (1 - e[inds]**2)

        return b


def time_since_periapsis(M, n, v, p, mu):
    if p.ndim == 0:
        if n>1e-19: # non-parabolic
            return M/n
        else: # parabolic
            D = tan(v/2)
            return sqrt(2*(p/2)**3/mu)*(D + D**3/3)
    else:
        parabolic = n<1e-19
        t = divide(M, n, out=zeros_like(p), where=~parabolic)

        D = tan(v[parabolic]/2)
        t[parabolic] = sqrt(2*(p[parabolic]/2)**3/mu)*(D + D**3/3)

        return t


def true_anomaly(e_vec, pos_vec, vel_vec, n_vec):
    if pos_vec.ndim == 1:
        if length_of(e_vec) > 1e-15: # not circular
            angle = angle_between(e_vec, pos_vec)
            v = angle if dots(pos_vec, vel_vec) > 0 else -angle%tau

        elif length_of(n_vec) < 1e-15: # circular and equatorial
            angle = arccos(pos_vec[0]/length_of(pos_vec))
            v = angle if vel_vec[0] < 0 else -angle%tau

        else: # circular and not equatorial
            angle = angle_between(n_vec, pos_vec)
            v = angle if pos_vec[2] >= 0 else -angle%tau

        return v if length_of(e_vec)<(1-1e-15) else normpi(v)
    else:
        v = zeros_like(pos_vec[0])
        circular = (length_of(e_vec)<1e-15)
        equatorial = (length_of(n_vec)<1e-15)

        inds = ~circular
        angle = angle_between(e_vec[:, inds], pos_vec[:, inds])
        condition = (dots(pos_vec[:, inds], vel_vec[:, inds]) > 0)
        v[inds] = where(condition, angle, -angle%tau)

        inds = circular*equatorial
        angle = arccos(pos_vec[0][inds]/length_of(pos_vec)[inds])
        condition = vel_vec[0][inds] < 0
        v[inds] = where(condition, angle, -angle%tau)

        inds = circular*~equatorial
        angle = angle_between(n_vec[:, inds], pos_vec[:, inds])
        condition = pos_vec[2][inds] >= 0
        v[inds] = where(condition, angle, -angle%tau)

        inds = length_of(e_vec)>(1-1e-15)
        v[inds] = normpi(v[inds])

        return v
