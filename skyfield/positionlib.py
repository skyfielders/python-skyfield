"""Classes representing different kinds of astronomical position."""

from numpy import array, einsum, exp

from .constants import RAD2DEG, TAU
from .data.spice import inertial_frames
#from .framelib import ICRS_to_J2000
from .functions import from_polar, length_of, to_polar, rot_z
from .earthlib import compute_limb_angle, refract
from .relativity import add_aberration, add_deflection
from .timelib import Time
from .units import Distance, Velocity, Angle, _interpret_angle

_ECLIPJ2000 = inertial_frames['ECLIPJ2000']
_GALACTIC = inertial_frames['GALACTIC']

class ICRS(object):
    """An x,y,z position whose axes are oriented to the ICRS system.

    The ICRS is a permanent coordinate system that has superseded the
    old series of equinox-based systems like B1900, B1950, and J2000.

    """
    geocentric = True
    altaz_rotation = None

    def __init__(self, position_au, velocity_au_per_d=None, jd=None):
        self.jd = jd
        self.position = Distance(position_au)
        if velocity_au_per_d is None:
            self.velocity = None
        else:
            self.velocity = Velocity(velocity_au_per_d)

    def __repr__(self):
        return '<%s position x,y,z au%s%s>' % (
            self.__class__.__name__,
            '' if (self.velocity is None) else
            ' and velocity xdot,ydot,zdot au/day',
            '' if self.jd is None else ' at date jd',
            )

    def __sub__(self, body):
        """Subtract two ICRS vectors to produce a third."""
        p = self.position.au - body.position.au
        if self.velocity is None or body.velocity is None:
            v = None
        else:
            v = body.velocity.au_per_d - self.velocity.au_per_d
        return ICRS(p, v, self.jd)

    def distance(self):
        """Return the length of this vector.

        >>> v = ICRS([1.0, 1.0, 0.0])
        >>> print(v.distance())
        1.41421 au

        """
        return Distance(length_of(self.position.au))

    def radec(self, epoch=None):
        """Return this position as a tuple (RA, declination, distance).

        >>> ra, dec, distance = ICRS([1.0, 1.0, 1.0]).radec()
        >>> ra
        <Angle 03h 00m 00.00s>
        >>> dec
        <Angle +35deg 15' 51.8">
        >>> distance
        <Distance 1.73205 au>

        """
        position_au = self.position.au
        if epoch is not None:
            if isinstance(epoch, Time):
                pass
            elif isinstance(epoch, float):
                epoch = Time(tt=epoch)
            elif epoch == 'date':
                epoch = self.jd
            else:
                raise ValueError('the epoch= must be a Julian date,'
                                 ' a floating point Terrestrial Time (TT),'
                                 ' or the string "date" for epoch-of-date')
            position_au = einsum('ij...,j...->i...', epoch.M, position_au)
        r_au, dec, ra = to_polar(position_au)
        return (Angle(radians=ra, preference='hours'),
                Angle(radians=dec, signed=True),
                Distance(r_au))

    def ecliptic_position(self):
        """Return an x,y,z position relative to the ecliptic plane."""
        vector = _ECLIPJ2000.dot(self.position.au)
        return Distance(vector)

    def ecliptic_latlon(self):
        """Return ecliptic latitude, longitude, and distance."""
        vector = _ECLIPJ2000.dot(self.position.au)
        d, lat, lon = to_polar(vector)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def galactic_position(self):
        """Return an x,y,z position relative to the galactic plane."""
        vector = _GALACTIC.dot(self.position.au)
        return Distance(vector)

    def galactic_latlon(self):
        """Return galactic latitude, longitude, and distance."""
        vector = _GALACTIC.dot(self.position.au)
        d, lat, lon = to_polar(vector)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def to_spice_frame(self, name):
        vector = self.position.au
        vector = inertial_frames[name].dot(vector)
        d, dec, ra = to_polar(vector)
        return (Angle(radians=ra, preference='hours', signed=True),
                Angle(radians=dec),
                Distance(au=d))

    def from_altaz(self, alt=None, az=None, alt_degrees=None, az_degrees=None):
        """Generate an Apparent position from an altitude and azimuth.

        The altitude and azimuth can each be provided as an `Angle`
        object, or else as a number of degrees provided as either a
        float or a tuple of degrees, arcminutes, and arcseconds::

            alt=Angle(...), az=Angle(...)
            alt_degrees=23.2289, az_degrees=142.1161
            alt_degrees=(23, 13, 44.1), az_degrees=(142, 6, 58.1)

        """
        # TODO: should this method live on another class?
        R = self.altaz_rotation
        if R is None:
            raise ValueError('only a position generated by a topos() call'
                             ' knows the orientation of the horizon'
                             ' and can understand altitude and azimuth')
        alt = _interpret_angle('alt', alt, alt_degrees)
        az = _interpret_angle('az', az, az_degrees)
        r = 0.1  # close enough to make gravitational refraction irrelevant
        p = from_polar(r, alt, az)
        p = einsum('ji...,j...->i...', R, p)
        return Apparent(p)


class Barycentric(ICRS):
    """BCRS: an ICRS x,y,z position measured from the Solar System barycenter.

    """
    def observe(self, body):
        """Return the astrometric position of `body` viewed from this position.

        """
        astrometric = body._observe_from_bcrs(self)
        return astrometric


class Astrometric(ICRS):
    """An astrometric position as an x,y,z vector in the ICRS.

    The *astrometric position* of a body is its position relative to an
    observer, adjusted for light-time delay: the position of the body
    back when it emitted (or reflected) the light that is now reaching
    the observer's eyes or telescope.  This is always a difference
    between two BCRS vectors.

    """
    def apparent(self):
        """Return the apparent position where this will appear in the sky.

        This method determines how relativity affects an image, and
        returns the :class:`~skyfield.positionlib.Apparent` position
        where the body will actually appear in the sky.  The effects
        modeled are the deflection that the image will experience if its
        light passes close to large masses in the Solar System, and the
        aberration caused by the observer's own velocity.

        These transforms convert the position from the BCRS reference
        frame of the Solar System barycenter and to the reference frame
        of the observer.  In the specific case of an Earth observer, the
        output reference frame is the GCRS.

        """
        jd = self.jd
        position_au = self.position.au.copy()
        observer = self.observer

        if observer.geocentric:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position_au, observer.position.au)
            include_earth_deflection = limb_angle >= 0.8

        add_deflection(position_au, observer.position.au, observer.ephemeris,
                       jd, include_earth_deflection)
        add_aberration(position_au, observer.velocity.au_per_d, self.light_time)

        a = Apparent(position_au, jd=jd)
        a.observer = self.observer
        return a


class Apparent(ICRS):
    """An apparent position as an x,y,z vector in the GCRS.

    The *apparent position* of a body is its position relative to an
    observer, adjusted not only for the light-time delay between the
    body and an observer (which was already accounted for in the
    object's astrometric position), but also adjusted for deflection
    (light rays bending as they pass large masses like the Sun or
    Jupiter) and aberration (light slanting because of the observer's
    motion through space).

    Included in aberration is the relativistic transformation that takes
    the position out of the BCRS centered on the solar system barycenter
    and into the GCRS centered on the Earth.

    If the observer was a planet or satellite with its own orbit around
    the Sun, then this apparent position is not really a GCRS position,
    but belongs to a GCRS-like system centered on that observer instead.

    """
    def altaz(self, temperature_C=None, pressure_mbar='standard'):
        """Return the position as a tuple ``(alt, az, distance)``.

        `alt` - Altitude in degrees above the horizon.
        `az` - Azimuth angle east around the horizon from due-north.
        `distance` - Distance to the object.

        """
        try:
            topos = self.observer.topos
            R = self.observer.altaz_rotation
        except AttributeError:
            raise ValueError('to compute an apparent position, you must'
                             ' observe from a specific Earth location that'
                             ' you specify using a Topos instance')

        # TODO: wobble

        position_au = einsum('ij...,j...->i...', R, self.position.au)
        r_au, alt, az = to_polar(position_au)

        if temperature_C is None:
            alt = Angle(radians=alt)
        else:
            if temperature_C == 'standard':
                temperature_C = 10.0
            if pressure_mbar == 'standard':
                pressure_mbar = 1010.0 * exp(-topos.elevation.m / 9.1e3)
            alt = refract(alt * RAD2DEG, temperature_C, pressure_mbar)
            alt = Angle(degrees=alt)

        return alt, Angle(radians=az), Distance(r_au)


class Geocentric(ICRS):
    """A position referred to the GCRS as measured from the geocenter."""

    def observe(self, other):
        gcrs_method = getattr(other, 'gcrs')
        if gcrs_method is None:
            raise ValueError('currently a Geocentric location can only'
                             ' observe an object that can generate a'
                             ' GCRS position through a .gcrs() method')
        g = gcrs_method(self.jd)
        # TODO: light-travel-time backdating, if distant enough?
        p = g.position.au - self.position.au
        v = g.velocity.au_per_d - self.velocity.au_per_d
        a = Apparent(p, v, self.jd)
        a.observer = self
        return a


def ITRF_to_GCRS(jd, rITRF):  # todo: velocity

    # Todo: wobble

    spin = rot_z(jd.gast * TAU / 24.0)
    position = einsum('ij...,j...->i...', spin, array(rITRF))
    return einsum('ij...,j...->i...', jd.MT, position)
