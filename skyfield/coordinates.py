"""Coordinate systems."""

import numpy as np
from numpy import arcsin, arctan2, ndarray, max, min, sqrt
from math import cos, sin, pi
from skyfield.angles import interpret_longitude, interpret_latitude
from skyfield.framelib import ICRS_to_J2000
from skyfield.nutationlib import compute_nutation
from skyfield.precessionlib import compute_precession
from skyfield.relativity import add_aberration, add_deflection

J2000 = 2451545.0
C_AUDAY = 173.1446326846693

ecliptic_obliquity = (23 + (26/60.) + (21.406/3600.)) * pi / 180.

class Hours(float):

    def hms(self):
        n = self if self > 0. else -self
        d, fraction = divmod(n * 12. / pi, 1.)
        m, s = divmod(fraction * 3600., 60.)
        return d if self > 0. else -d, m, s

    def __str__(self):
        return '%d:%02d:%02f' % self.hms()

class Degrees(float):

    def dms(self):
        n = self if self > 0. else -self
        d, fraction = divmod(n * 180. / pi, 1.)
        m, s = divmod(fraction * 3600., 60.)
        return d if self > 0. else -d, m, s

    def __str__(self):
        return '%d:%02d:%02f' % self.dms()

class XYZ(object):

    def __init__(self, position, velocity=None, jd=None):
        self.position = position
        self.velocity = velocity
        self.jd = jd

class ICRS(XYZ):

    geocentric = True

    def observe(self, body):
        # TODO: should also accept another ICRS?

        jd = self.jd
        lighttime0 = 0.0
        vector = body(jd).position - self.position
        euclidian_distance = distance = sqrt((vector * vector).sum(axis=0))

        for i in range(10):
            lighttime = distance / C_AUDAY
            delta = lighttime - lighttime0
            if -1e-12 < min(delta) and max(delta) < 1e-12:
                break
            lighttime0 = lighttime
            target = body(jd - lighttime)
            vector = target.position - self.position
            distance = sqrt((vector * vector).sum(axis=0))
        else:
            raise ValueError('observe() light-travel time failed to converge')

        g = GCRS(vector, target.velocity - self.velocity, jd)
        g.observer = self
        g.distance = euclidian_distance
        g.lighttime = lighttime
        return g

class Topos(object):

    def __init__(self, longitude, latitude, elevation=0.,
                 temperature=10.0, pressure=1010.0):
        self.longitude = interpret_longitude(longitude)
        self.latitude = interpret_latitude(latitude)
        self.elevation = elevation

    def __call__(self, jd_tt):
        from skyfield.earthlib import geocentric_position_and_velocity
        if not hasattr(jd_tt, 'shape'):
            jd_tt = np.array((jd_tt,))
        e = self.earth(jd_tt)
        tpos, tvel = geocentric_position_and_velocity(self, jd_tt)
        t = ToposICRS(e.position + tpos, e.velocity + tvel, jd_tt)
        t.ephemeris = self.ephemeris
        return t

class ToposICRS(ICRS):
    """In ICRS, right?"""

    geocentric = False

class GCRS(XYZ):

    def astrometric(self, epoch=None):
        eq = Astrometric()
        eq.ra, eq.dec = to_polar(self.position)
        eq.distance = self.distance
        # TODO: epoch
        return eq

    def apparent(self):
        jd = self.jd
        position = self.position.copy()

        observer = self.observer

        from skyfield.earthlib import compute_limb_angle

        if observer.geocentric:
            use_earth = np.array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position, observer.position)
            use_earth = limb_angle >= 0.8

        add_deflection(position, observer.position, observer.ephemeris,
                       jd, use_earth)
        add_aberration(position, observer.velocity, self.lighttime)

        if np.isscalar(jd):
            jd = np.array((jd,))
            position = position.reshape((3, 1))
        else:
            position = position.copy()

        position = position.T.dot(ICRS_to_J2000)
        position = np.einsum('ij,jki->ik', position, compute_precession(jd))
        position = np.einsum('ij,jki->ik', position, compute_nutation(jd))
        position = position.T

        eq = Apparent()
        eq.ra, eq.dec = to_polar(position)
        eq.distance = self.distance
        return eq

class RADec():
    def __repr__(self):
        return '<%s position RA=%s dec=%s>' % (
            self.__class__.__name__, self.ra, self.dec)

class Astrometric(RADec):
    pass

class Apparent(RADec):
    pass

class HeliocentricLonLat(ndarray):

    def __new__(cls, other):
        self = ndarray.__new__(cls, (3,))
        if isinstance(other, ICRS):
            x, y, z = other
            y, z = (
                y * cos(ecliptic_obliquity) + z * sin(ecliptic_obliquity),
                z * cos(ecliptic_obliquity) - y * sin(ecliptic_obliquity),
                )
            self[2] = r = sqrt(x*x + y*y + z*z)
            self[0] = arctan2(-y, -x) + pi
            self[1] = arcsin(z / r)
        else:
            raise ValueError('how do I use that?')
        return self

    @property
    def lon(self): return Degrees(self[0])

    @property
    def lat(self): return Degrees(self[1])

    @property
    def r(self): return self[2]

def to_polar(position):
    r = sqrt((position * position).sum(axis=0))
    phi = arctan2(-position[1], -position[0]) + pi
    theta = arcsin(position[2] / r)
    return phi, theta
