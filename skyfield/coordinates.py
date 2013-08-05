"""Coordinate systems."""

from numpy import (arcsin, arctan2, array, cos, einsum, isscalar,
                   ndarray, pi, sin, sqrt)
from skyfield.angles import interpret_longitude, interpret_latitude, Angle, \
    HourAngle, tau
from skyfield.framelib import ICRS_to_J2000
from skyfield.nutationlib import compute_nutation
from skyfield.precessionlib import compute_precession
from skyfield.relativity import add_aberration, add_deflection

J2000 = 2451545.0
C_AUDAY = 173.1446326846693

ecliptic_obliquity = (23 + (26/60.) + (21.406/3600.)) * pi / 180.


class XYZ(object):

    def __init__(self, position, velocity=None, jd=None):
        self.jd = jd
        if getattr(jd, 'shape', ()) == ():
            self.position = position.reshape((3,))
            self.velocity = velocity.reshape((3,))
        else:
            self.position = position
            self.velocity = velocity

    def __repr__(self):
        return '<%s position x,y,z AU%s%s>' % (
            self.__class__.__name__,
            '' if (self.velocity is None) else
            ' and velocity xdot,ydot,zdot AU/day',
            '' if self.jd is None else ' at date jd',
            )

    @property
    def x(self): return self.position[0]

    @property
    def y(self): return self.position[1]

    @property
    def z(self): return self.position[2]

    @property
    def xdot(self): return self.velocity[0]

    @property
    def ydot(self): return self.velocity[1]

    @property
    def zdot(self): return self.velocity[2]

class ICRS(XYZ):

    geocentric = True

    def __sub__(self, body):
        if self.velocity is None or body.velocity is None:
            velocity = None
        else:
            velocity = body.velocity - self.velocity
        return XYZ(self.position - body.position, velocity, self.jd)

class Topos(object):

    def __init__(self, longitude, latitude, elevation=0.,
                 temperature=10.0, pressure=1010.0):
        self.longitude = interpret_longitude(longitude)
        self.latitude = interpret_latitude(latitude)
        self.elevation = elevation

    def __call__(self, jd_tt):
        from skyfield.earthlib import geocentric_position_and_velocity
        if not hasattr(jd_tt, 'shape'):
            jd_tt = array((jd_tt,))
        e = self.earth(jd_tt)
        tpos, tvel = geocentric_position_and_velocity(self, jd_tt)
        t = ToposICRS(e.position + tpos, e.velocity + tvel, jd_tt)
        t.latitude = self.latitude
        t.longitude = self.longitude
        t.ephemeris = self.ephemeris
        return t

class ToposICRS(ICRS):
    """In ICRS, right?"""

    geocentric = False

class GCRS(XYZ):

    def astrometric(self, epoch=None):
        eq = Astrometric()
        eq.ra, eq.dec = to_polar(self.position, phi_class=HourAngle)
        eq.position = self.position
        eq.distance = self.distance
        eq.observer = self.observer
        # TODO: epoch
        return eq

    def apparent(self):
        jd = self.jd
        position = self.position.copy()

        observer = self.observer

        from skyfield.earthlib import compute_limb_angle

        if observer.geocentric:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position, observer.position)
            include_earth_deflection = limb_angle >= 0.8

        add_deflection(position, observer.position, observer.ephemeris,
                       jd, include_earth_deflection)
        add_aberration(position, observer.velocity, self.lighttime)

        if isscalar(jd):
            jd = array((jd,))
            position = position.reshape((3, 1))
        else:
            position = position.copy()

        position = position.T.dot(ICRS_to_J2000)
        position = einsum('ij,jki->ik', position, compute_precession(jd))
        position = einsum('ij,jki->ik', position, compute_nutation(jd))
        position = position.T

        eq = Apparent()
        eq.ra, eq.dec = to_polar(position, phi_class=HourAngle)
        eq.jd = jd
        eq.position = position
        eq.distance = self.distance
        eq.observer = self.observer
        return eq

class RADec(object):
    def __repr__(self):
        return '<%s position RA=%r dec=%r>' % (
            self.__class__.__name__, self.ra, self.dec)

class Astrometric(RADec):
    pass

class Apparent(RADec):
    """Topocentric RA and declination vs true equator and equinox-of-date."""

    def horizontal(self):
        lat = self.observer.latitude
        lon = self.observer.longitude
        delta_t = 0.0

        sinlat = sin(lat)
        coslat = cos(lat)
        sinlon = sin(lon)
        coslon = cos(lon)

        uze = array([coslat * coslon, coslat * sinlon, sinlat])
        une = array([-sinlat * coslon, -sinlat * sinlon, coslat])
        uwe = array([sinlon, -coslon, 0.0])

        # TODO: allow called to ask for corrections using xp and yp

        # TODO: this might be indefensible if delta_t is not zero
        jd_ut1 = self.jd

        from .timescales import sidereal_time

        def spin(angle, vector):
            cosang = cos(angle)
            sinang = sin(angle)

            r1 = array([cosang, sinang, 0.0])
            r2 = array([-sinang, cosang, 0.0])
            r3 = array([0.0, 0.0, 1.0])

            return array([
                r1.dot(vector),
                r2.dot(vector),
                r3.dot(vector),
                ])

        gast = sidereal_time(jd_ut1, delta_t, use_eqeq=True)
        uz = spin(-gast * tau / 24.0, uze)
        un = spin(-gast * tau / 24.0, une)
        uw = spin(-gast * tau / 24.0, uwe)

        p = self.position[:,0]  # TODO: vectorize this whole operation
        pz = p.dot(uz)
        pn = p.dot(un)
        pw = p.dot(uw)

        proj = sqrt(pn * pn + pw * pw)

        h = Horizontal()

        if proj > 0.0:
            h.az = -arctan2(pw, pn)

        if h.az < 0.0:
            h.az += tau

        if h.az >= tau:
            h.az -= tau

        h.zd = arctan2(proj, pz)
        h.distance = self.distance
        # TODO: which JD to save with coordinate?
        return h

class Horizontal(object):
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

    # @property
    # def lon(self): return Degrees(self[0])

    # @property
    # def lat(self): return Degrees(self[1])

    @property
    def r(self): return self[2]

def to_polar(position, phi_class=Angle):
    r = sqrt((position * position).sum(axis=0))
    phi = phi_class(r.shape)
    theta = Angle(r.shape)
    arctan2(-position[1], -position[0], out=phi)
    arcsin(position[2] / r, out=theta)
    phi += pi
    return phi, theta
