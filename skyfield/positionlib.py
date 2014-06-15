"""Classes representing different kinds of astronomical position."""

from numpy import arcsin, arctan2, array, cos, einsum, pi, sin

from .constants import TAU
from .functions import length_of, spin_x
from .earthlib import (compute_limb_angle, geocentric_position_and_velocity,
                       sidereal_time)
from .functions import dots
from .relativity import add_aberration, add_deflection
from .timelib import JulianDate, takes_julian_date
from .units import (Distance, Velocity, Angle, HourAngle, SignedAngle,
                    _interpret_ltude)

ecliptic_obliquity = (23 + (26/60.) + (21.406/3600.)) * pi / 180.
quarter_tau = 0.25 * TAU

class ICRS(object):
    """An x,y,z position whose axes are oriented to the ICRS system.

    The ICRS is a permanent coordinate system that has superseded the
    old series of equinox-based systems like B1900, B1950, and J2000.

    """
    geocentric = True  # TODO: figure out what this meant and get rid of it

    def __init__(self, position_AU, velocity_AU_per_d=None, jd=None):
        self.jd = jd
        self.position = Distance(position_AU)
        if velocity_AU_per_d is None:
            self.velocity = None
        else:
            self.velocity = Velocity(velocity_AU_per_d)

    def __repr__(self):
        return '<%s position x,y,z AU%s%s>' % (
            self.__class__.__name__,
            '' if (self.velocity is None) else
            ' and velocity xdot,ydot,zdot AU/day',
            '' if self.jd is None else ' at date jd',
            )

    def __sub__(self, body):
        """Subtract two ICRS vectors to produce a third."""
        p = self.position.AU - body.position.AU
        if self.velocity is None or body.velocity is None:
            v = None
        else:
            v = body.velocity.AU_per_d - self.velocity.AU_per_d
        return ICRS(p, v, self.jd)

    def observe(self, body):
        return body.observe_from(self)

    def radec(self, epoch=None):
        position_AU = self.position.AU
        if epoch is not None:
            if isinstance(epoch, JulianDate):
                pass
            elif isinstance(epoch, float):
                epoch = JulianDate(tt=epoch)
            elif epoch == 'date':
                epoch = self.jd
            else:
                raise ValueError('the epoch= must be a Julian date,'
                                 ' a floating point Terrestrial Time (TT),'
                                 ' or the string "date" for epoch-of-date')
            position_AU = einsum('ij...,j...->i...', epoch.M, position_AU)
        r_AU, dec, ra = to_polar(position_AU)
        return HourAngle(radians=ra), SignedAngle(radians=dec), Distance(r_AU)

# class to represent a point in the IC reference frame
class ICRCoordinates:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def equalTo(self, other):
        # TODO: override ==, and add epsilons here
        return (self.x == other.x) and (self.y == other.y) and (self.z == other.z)

    def __repr__(self):
        return '(%s, %s, %s)' % (self.x, self.y, self.z)

class Topos(object):

    def __init__(self, latitude=None, longitude=None,
                 latitude_degrees=None, longitude_degrees=None,
                 elevation=0., temperature=10.0, pressure=1010.0):

        if latitude_degrees is not None:
            latitude = Angle(degrees=latitude_degrees)
        elif isinstance(latitude, str):
            latitude = _interpret_ltude(latitude, 'latitude', 'N', 'S')
        elif not isinstance(latitude, Angle):
            raise TypeError('please provide either latitude_degrees=<float>'
                            ' or else latitude=<skyfield.units.Angle object>'
                            ' with north being positive')

        if longitude_degrees is not None:
            longitude = Angle(degrees=longitude_degrees)
        elif isinstance(longitude, str):
            longitude = _interpret_ltude(longitude, 'longitude', 'E', 'W')
        elif not isinstance(longitude, Angle):
            raise TypeError('please provide either longitude_degrees=<float>'
                            ' or else longitude=<skyfield.units.Angle object>'
                            ' with east being positive')

        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation

        lat = latitude.radians()
        lon = longitude.radians()
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
        tpos_AU, tvel_AU_per_d = geocentric_position_and_velocity(self, jd)
        t = ToposICRS(e.position.AU + tpos_AU,
                      e.velocity.AU_per_d + tvel_AU_per_d,
                      jd)
        t.rGCRS = tpos_AU
        t.vGCRS = tvel_AU_per_d
        t.topos = self
        t.ephemeris = self.ephemeris
        return t

class ToposICRS(ICRS):
    """In ICRS, right?"""

    geocentric = False

class Barycentric(ICRS):
    """An ICRS x,y,z position referenced to the Solar System barycenter."""

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
        position_AU = self.position.AU.copy()
        observer = self.observer

        if observer.geocentric:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position_AU, observer.position.AU)
            include_earth_deflection = limb_angle >= 0.8

        add_deflection(position_AU, observer.position.AU, observer.ephemeris,
                       jd.tdb, include_earth_deflection)
        add_aberration(position_AU, observer.velocity.AU_per_d, self.lighttime)

        a = Apparent(position_AU, jd=jd)
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

        p = einsum('ij...,j...->i...', self.jd.M, self.position.AU)

        pz = dots(p, uz)
        pn = dots(p, un)
        pw = dots(p, uw)

        position_AU = array([pn, -pw, pz])

        r_AU, alt, az = to_polar(position_AU)
        return Angle(radians=alt), Angle(radians=az), Distance(r_AU)


def to_polar(xyz):
    """Convert ``[x y z]`` into spherical coordinates ``(r, theta, phi)``.

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    The order of the three return values is intended to match ISO 31-11.

    """
    r = length_of(xyz)
    theta = arcsin(xyz[2] / r)
    phi = arctan2(-xyz[1], -xyz[0]) + pi
    return r, theta, phi

def ITRF_to_GCRS(jd, rITRF):  # todo: velocity

    # Todo: wobble

    gast = sidereal_time(jd, use_eqeq=True)
    spin = spin_x(-gast * TAU / 24.0)
    position = einsum('i...,ij...->j...', array(rITRF), spin)
    return einsum('ij...,j...->i...', jd.MT, position)
