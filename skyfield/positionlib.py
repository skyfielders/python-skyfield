"""Classes representing different kinds of astronomical position."""

from numpy import arcsin, arctan2, array, cos, einsum, ndarray, pi, sin, sqrt

from .angles import interpret_longitude, interpret_latitude, Angle, HourAngle
from .constants import TAU
from .functions import length_of, spin_x
from .earthlib import (compute_limb_angle, geocentric_position_and_velocity,
                       sidereal_time)
from .framelib import ICRS_to_J2000, J2000_to_ICRS
from .functions import dots
from .relativity import add_aberration, add_deflection
from .timescales import takes_julian_date

ecliptic_obliquity = (23 + (26/60.) + (21.406/3600.)) * pi / 180.
quarter_tau = 0.25 * TAU

class ICRS(object):
    """An x,y,z position whose axes are oriented to the ICRS system.

    The ICRS is a permanent coordinate system that has superseded the
    old series of equinox-based systems like B1900, B1950, and J2000.

    """
    geocentric = True  # TODO: figure out what this meant and get rid of it

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
        d, dec, ra = to_polar(position, phi_class=HourAngle)
        return ra, dec, d

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
        t.rGCRS = tpos
        t.vGCRS = tvel
        t.topos = self
        t.ephemeris = self.ephemeris
        return t

class ToposICRS(ICRS):
    """In ICRS, right?"""

    geocentric = False

class Astrometric(ICRS):
    """An astrometric position as an x,y,z vector in the ICRS.

    The *astrometric position* of a body is its position relative to an
    observer, adjusted for light-time delay: the position of the body
    back when it emitted (or reflected) the light or other radiation
    that is just now reaching the observer's eyes or telescope.  This is
    always a difference between two BCRS vectors.

    """
    def apparent(self):
        """Return the corresponding apparent position."""
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

        a = Apparent(position, jd=jd)
        a.observer = self.observer
        return a

class Apparent(ICRS):
    """An apparent position as an x,y,z vector in the GCRS.

    The *apparent position* of a body is its position relative to an
    observer, adjusted not only for the light-time delay between the
    body and an observer (which was already accounted for in the
    object's astrometric position), but also adjusted for deflection
    (its light rays bending as they pass large masses like the Sun or
    Jupiter) and aberration (light slanting because of the observer's
    motion through space).

    Included in aberration is the relativistic transformation that takes
    the position out of the BCRS centered on the solar system barycenter
    and into the GCRS centered on the Earth.

    If the observer was a planet or satellite with its own orbit around
    the Sun, then this apparent position is not really a GCRS position,
    but belongs to a GCRS-like system centered on that observer instead.

    """
    def altaz(self):
        """Return the position as a tuple ``(alt, az, d)``.

        `alt` - Altitude in degrees above the horizon.
        `az` - Azimuth angle east around the horizon from due-north.
        `d` - Distance to the object.

        """
        try:
            topos = self.observer.topos
            uze = topos.up
            une = topos.north
            uwe = topos.west
        except AttributeError:
            raise ValueError('to compute an apparent position, you must'
                             ' observe from a specific Earth location that'
                             ' you specify using a Topos instance')

        # TODO: wobble

        gast = sidereal_time(self.jd, use_eqeq=True)
        spin = spin_x(-gast * TAU / 24.0)
        uz = einsum('i...,ij...->j...', uze, spin)
        un = einsum('i...,ij...->j...', une, spin)
        uw = einsum('i...,ij...->j...', uwe, spin)

        p = self.position
        p = to_epoch(p, self.jd)

        pz = dots(p, uz)
        pn = dots(p, un)
        pw = dots(p, uw)

        position = array([pn, -pw, pz])

        d, alt, az = to_polar(position)
        return alt, az, d


def to_polar(position, phi_class=Angle):
    """Convert ``[x y z]`` into spherical coordinates ``(r, theta, phi)``.

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    The order of values in the tuple is intended to match ISO 31-11.

    """
    r = length_of(position)
    theta = Angle(radians=arcsin(position[2] / r))
    phi = phi_class(radians=arctan2(-position[1], -position[0]) + pi)
    return r, theta, phi

def to_epoch(position, epoch):
    """Convert an ICRS position into an Earth equatorial position."""

    jd = epoch

    position = position.T.dot(ICRS_to_J2000)
    position = einsum('...j,jk...->...k', position, jd.PT)
    position = einsum('...j,jk...->...k', position, jd.NT)
    position = position.T

    return position

def ITRF_to_GCRS(jd, rITRF):  # todo: velocity

    # Todo: wobble

    gast = sidereal_time(jd, use_eqeq=True)
    spin = spin_x(-gast * TAU / 24.0)
    position = einsum('i...,ij...->j...', array(rITRF), spin)

    position = array(position)
    position = position.T

    position = einsum('...j,jk...->...k', position, jd.N)
    position = einsum('...j,jk...->...k', position, jd.P)
    position = position.dot(J2000_to_ICRS)
    rGCRS = position.T
    return rGCRS
