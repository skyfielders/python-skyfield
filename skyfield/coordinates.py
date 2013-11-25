"""Coordinate systems."""

from numpy import (arcsin, arctan2, array, cos, einsum,
                   ndarray, ones_like, pi, sin, sqrt, zeros_like)

from .angles import (interpret_longitude, interpret_latitude, Angle, HourAngle)
from .constants import TAU
from .functions import length
from .earthlib import (compute_limb_angle, geocentric_position_and_velocity,
                       sidereal_time)
from .framelib import ICRS_to_J2000
from .functions import dots
from .nutationlib import compute_nutation
from .precessionlib import compute_precession
from .relativity import add_aberration, add_deflection
from .timescales import takes_julian_date

ecliptic_obliquity = (23 + (26/60.) + (21.406/3600.)) * pi / 180.
quarter_tau = 0.25 * TAU

class XYZ(object):

    def __init__(self, position, velocity=None, jd=None):
        self.jd = jd
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
    """An x,y,z position whose axes are oriented to the ICRS system.

    The ICRS is a permanent coordinate system that has superseded the
    old series of equinox-based systems like B1900, B1950, and J2000.

    """
    geocentric = True  # TODO: figure out what this meant and get rid of it

    def __sub__(self, body):
        """Subtract two ICRS vectors to produce a third."""
        if self.velocity is None or body.velocity is None:
            velocity = None
        else:
            velocity = body.velocity - self.velocity
        return ICRS(self.position - body.position, velocity, self.jd)

    def observe(self, body):
        return body.observe_from(self)

    def radec(self, epoch=None):
        position = self.position
        if epoch is not None:
            position = to_epoch(position, self.jd)
        return to_polar(position, phi_class=HourAngle)

class Topos(object):

    def __init__(self, longitude, latitude, elevation=0.,
                 temperature=10.0, pressure=1010.0):
        self.longitude = lon = interpret_longitude(longitude)
        self.latitude = lat = interpret_latitude(latitude)
        self.elevation = elevation

        sinlat = sin(lat)
        coslat = cos(lat)
        sinlon = sin(lon)
        coslon = cos(lon)

        self.up = array([coslat * coslon, coslat * sinlon, sinlat])
        self.north = array([-sinlat * coslon, -sinlat * sinlon, coslat])
        self.west = array([sinlon, -coslon, 0.0])

    @takes_julian_date
    def __call__(self, jd):
        """Compute where this Earth location was in space on a given date."""
        e = self.ephemeris.earth(jd)
        tpos, tvel = geocentric_position_and_velocity(self, jd)
        t = ToposICRS(e.position + tpos, e.velocity + tvel, jd)
        t.topos = self
        t.ephemeris = self.ephemeris
        return t

class ToposICRS(ICRS):
    """In ICRS, right?"""

    geocentric = False

class Astrometric(ICRS):
    """An astrometric position as GCRS x,y,z coordinates.

    The *astrometric position* of a body is its position adjusted for
    the light-time delay between the body and an observer; it is the
    position of the body back when it emitted (or reflected) the light
    that is just now reaching the observer's eyes or telescope.

    """
    def apparent(self):
        jd = self.jd
        position = self.position.copy()
        observer = self.observer

        if observer.geocentric:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position, observer.position)
            include_earth_deflection = limb_angle >= 0.8

        add_deflection(position, observer.position, observer.ephemeris,
                       jd.tdb, include_earth_deflection)
        add_aberration(position, observer.velocity, self.lighttime)

        # p = compute_precession(jd_tdb)
        # n = compute_nutation(jd)

        # position = position.T.dot(ICRS_to_J2000)
        # position = einsum('...j,jk...->...k', position, p)
        # position = einsum('...j,jk...->...k', position, n)
        # position = position.T

        a = Apparent(position, jd=jd)
        a.distance = self.distance
        a.observer = self.observer
        return a

class Apparent(ICRS):
    """An apparent position as GCRS x,y,z coordinates.

    The *apparent position* of a body is its position adjusted not only
    for the light-time delay between the body and an observer, which is
    already accounted for in the object's astrometric position, but
    further adjusted for aberration (Earth's motion through slants the
    light rays that we see from its the surface) and deflection (light
    bends as it passes large masses like the Sun or Jupiter).

    """
    def horizontal(self):

        try:
            topos = self.observer.topos
            uze = topos.up
            une = topos.north
            uwe = topos.west
        except AttributeError:
            raise ValueError('to compute an apparent position, you must'
                             ' observe from a specific Earth location that'
                             ' you specify using a Topos instance')

        # TODO: allow for corrections using xp and yp

        def spin(angle, vector):
            z = zeros_like(angle)
            u = ones_like(angle)
            cosang = cos(angle)
            sinang = sin(angle)

            rotation = array([
                [cosang, -sinang, z],
                [sinang, cosang, z],
                [z, z, u],
                ])

            return einsum('i...,ij...->j...', vector, rotation)

        gast = sidereal_time(self.jd, use_eqeq=True)
        uz = spin(-gast * TAU / 24.0, uze)
        un = spin(-gast * TAU / 24.0, une)
        uw = spin(-gast * TAU / 24.0, uwe)

        p = self.position
        p = to_epoch(p, self.jd)

        pz = dots(p, uz)
        pn = dots(p, un)
        pw = dots(p, uw)

        position = array([pn, -pw, pz])

        h = Horizontal()
        h.az, h.alt, d = to_polar(position)
        h.zd = quarter_tau - h.alt
        h.jd = self.jd
        h.distance = self.distance
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
    r = length(position)
    phi = phi_class(r.shape)
    theta = Angle(r.shape)
    arctan2(-position[1], -position[0], out=phi)
    arcsin(position[2] / r, out=theta)
    phi += pi
    return phi, theta, r

def to_epoch(position, epoch):
    jd = epoch
    p = compute_precession(jd.tdb)
    n = compute_nutation(jd)

    position = position.T.dot(ICRS_to_J2000)
    position = einsum('...j,jk...->...k', position, p)
    position = einsum('...j,jk...->...k', position, n)
    position = position.T

    return position
