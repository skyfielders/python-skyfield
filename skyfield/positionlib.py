"""Classes representing different kinds of astronomical position."""

from numpy import array, arccos, clip, einsum, exp

from .constants import RAD2DEG, tau
from .data.spice import inertial_frames
from .functions import dots, from_polar, length_of, to_polar, rot_z
from .earthlib import compute_limb_angle, refract
from .relativity import add_aberration, add_deflection
from .timelib import Time
from .units import Distance, Velocity, Angle, _interpret_angle

_ECLIPJ2000 = inertial_frames['ECLIPJ2000']
_GALACTIC = inertial_frames['GALACTIC']


def build_position(position_au, velocity_au_per_d=None, t=None,
                   center=None, target=None, observer_data=None):
    if center == 0:
        cls = Barycentric
    elif center == 399:
        cls = Geocentric
    elif observer_data is not None:
        # TODO: is the presence of observer_data enough to justify
        # Apparent?  Or could that in some cases skip the important
        # corrections that take place in the .apparent() method?
        cls = Apparent
    else:
        cls = ICRF
    return cls(position_au, velocity_au_per_d, t, center, target, observer_data)


class ICRF(object):
    """An (x, y, z) position and velocity oriented to the ICRF axes.

    The ICRF is a permanent coordinate system that has superseded the
    old series of equinox-based systems like B1900, B1950, and J2000.

    """
    altaz_rotation = None

    def __init__(self, position_au, velocity_au_per_d=None, t=None,
                 center=None, target=None, observer_data=None):
        self.t = t
        self.position = Distance(position_au)
        if velocity_au_per_d is None:
            self.velocity = None
        else:
            self.velocity = Velocity(velocity_au_per_d)
        self.center = center
        self.target = target
        self.observer_data = observer_data

    def __repr__(self):
        return '<{0} position{1}{2}{3}{4}>'.format(
            self.__class__.__name__,
            '' if (self.velocity is None) else ' and velocity',
            '' if self.t is None else ' at date t',
            '' if self.center is None else ' center={0}'.format(self.center),
            '' if self.target is None else ' target={0}'.format(self.target),
        )

    def __sub__(self, body):
        """Subtract two ICRF vectors to produce a third."""
        p = self.position.au - body.position.au
        if self.velocity is None or body.velocity is None:
            v = None
        else:
            v = body.velocity.au_per_d - self.velocity.au_per_d
        return ICRF(p, v, self.t)

    def distance(self):
        """Compute the distance from the origin to this position.

        >>> v = ICRF([1, 1, 0])
        >>> print(v.distance())
        1.41421 au

        """
        return Distance(length_of(self.position.au))

    def speed(self):
        """Compute the magnitude of the velocity vector.

        >>> v = ICRF([0, 0, 0], [1, 2, 3])
        >>> print(v.speed())
        3.74166 au/day

        """
        return Velocity(length_of(self.velocity.au_per_d))

    def radec(self, epoch=None):
        r"""Compute equatorial (RA, declination, distance)

        When called without a parameter, this returns standard ICRF
        right ascension and declination:

        >>> ra, dec, distance = ICRF([1, 2, 3]).radec()
        >>> print(ra, dec, distance, sep='\n')
        04h 13m 44.39s
        +53deg 18' 02.8"
        3.74166 au

        If you instead want the coordinates referenced to the dynamical
        system defined by the Earth's mean equator and equinox, provide
        an epoch time.  To get J2000.0 coordinates, for example:

        >>> ra, dec, distance = ICRF([1, 2, 3]).radec(ts.J2000)
        >>> print(ra, dec, sep='\n')
        04h 13m 43.32s
        +53deg 17' 55.1"

        """
        position_au = self.position.au
        if epoch is not None:
            if isinstance(epoch, Time):
                pass
            elif isinstance(epoch, float):
                epoch = Time(None, tt=epoch)
            elif epoch == 'date':
                epoch = self.t
            else:
                raise ValueError('the epoch= must be a Time object,'
                                 ' a floating point Terrestrial Time (TT),'
                                 ' or the string "date" for epoch-of-date')
            position_au = einsum('ij...,j...->i...', epoch.M, position_au)
        r_au, dec, ra = to_polar(position_au)
        return (Angle(radians=ra, preference='hours'),
                Angle(radians=dec, signed=True),
                Distance(r_au))

    def separation_from(self, another_icrf):
        """Return the angle between this position and another.

        >>> print(ICRF([1,0,0]).separation_from(ICRF([1,1,0])))
        45deg 00' 00.0"

        You can also compute separations across an array of positions.

        >>> directions = ICRF([[1,0,-1,0], [0,1,0,-1], [0,0,0,0]])
        >>> directions.separation_from(ICRF([0,1,0])).degrees
        array([  90.,    0.,   90.,  180.])

        """
        p1 = self.position.au
        p2 = another_icrf.position.au
        u1 = p1 / length_of(p1)
        u2 = p2 / length_of(p2)
        if u2.ndim > 1:
            if u1.ndim == 1:
                u1 = u1[:,None]
        elif u1.ndim > 1:
            u2 = u2[:,None]
        c = dots(u1, u2)
        return Angle(radians=arccos(clip(c, -1.0, 1.0)))

    def ecliptic_position(self):
        """Compute J2000 ecliptic coordinates (x, y, z)"""
        vector = _ECLIPJ2000.dot(self.position.au)
        return Distance(vector)

    def ecliptic_latlon(self):
        """Compute J2000 ecliptic coordinates (lat, lon, distance)"""
        vector = _ECLIPJ2000.dot(self.position.au)
        d, lat, lon = to_polar(vector)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def galactic_position(self):
        """Compute galactic coordinates (x, y, z)"""
        vector = _GALACTIC.dot(self.position.au)
        return Distance(vector)

    def galactic_latlon(self):
        """Compute galactic coordinates (lat, lon, distance)"""
        vector = _GALACTIC.dot(self.position.au)
        d, lat, lon = to_polar(vector)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def _to_spice_frame(self, name):
        """Return (ra, dec, """
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
        # TODO: remove this `if` statement
        if self.observer_data:
            R = self.observer_data.altaz_rotation
        else:
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


# For compatibility with my original name for the class.  Not an
# important enough change to warrant a deprecation error for users, so:
ICRS = ICRF


class Barycentric(ICRF):
    """An (x, y, z) position measured from the Solar System barycenter.

    Each barycentric position is an ICRS position vector, meaning that
    the coordinate axes are defined by the high-precision ICRF that has
    replaced the old J2000.0 reference frame, and the coordinate origin
    is the BCRS gravitational center of the Solar System.

    Skyfield generates a `Barycentric` position whenever you ask a Solar
    System body for its location at a particular time:

    >>> t = ts.utc(2003, 8, 29)
    >>> mars.at(t)
    <Barycentric position and velocity at date t>

    """
    # Positions on the Earth's surface or in low Earth orbit can set
    # these GCRS vectors if they want apparent positions to include the
    # effect of the Earth's gravitational deflection of light.
    # TODO: move these on to observer_info
    _gcrs_position = None
    _gcrs_velocity = None

    def observe(self, body):
        """Compute the `Astrometric` position of a body from this location.

        To compute the body's astrometric position, it is first asked
        for its position at the time `t` of this position itself.  The
        distance to the body is then divided by the speed of light to
        find how long it takes its light to arrive.  Finally, the light
        travel time is subtracted from `t` and the body is asked for a
        series of increasingly exact positions to learn where it was
        when it emitted the light that is now reaching this position.

        >>> earth.at(t).observe(mars)
        <Astrometric position and velocity at date t>

        """
        # TODO: can _observe_from_bcrs() just return some vectors?
        astrometric = body._observe_from_bcrs(self)
        astrometric.observer_data = self.observer_data
        return astrometric


class Astrometric(ICRF):
    """An astrometric (x, y, z) position relative to a particular observer.

    The *astrometric position* of a body is its position relative to an
    observer, adjusted for light-time delay: the position of the body
    back when it emitted (or reflected) the light that is now reaching
    the observer's eyes or telescope.

    Astrometric positions are usually generated in Skyfield by calling
    the `Barycentric` method `observe()` to determine where a body will
    appear in the sky relative to a specific observer.

    """
    def apparent(self):
        """Compute an :class:`Apparent` position for this body.

        This applies two effects to the position that arise from
        relativity and shift slightly where the other body will appear
        in the sky: the deflection that the image will experience if its
        light passes close to large masses in the Solar System, and the
        aberration of light caused by the observer's own velocity.

        >>> earth.at(t).observe(mars).apparent()
        <Apparent position at date t>

        These transforms convert the position from the BCRS reference
        frame of the Solar System barycenter and to the reference frame
        of the observer.  In the specific case of an Earth observer, the
        output reference frame is the GCRS.

        """
        t = self.t
        position_au = self.position.au.copy()
        # TODO: get rid of self.observer
        observer = self.observer

        gcrs_position = (observer._gcrs_position
                         or self.observer_data.gcrs_position)

        if gcrs_position is None:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position_au, gcrs_position)
            include_earth_deflection = nadir_angle >= 0.8

        add_deflection(position_au, observer.position.au,
                       self.observer_data.ephemeris, t,
                       include_earth_deflection)

        add_aberration(position_au, observer.velocity.au_per_d, self.light_time)

        apparent = Apparent(position_au, t=t)
        apparent.observer = self.observer
        apparent.observer_data = self.observer_data
        return apparent


class Apparent(ICRF):
    """An apparent (x, y, z) position relative to a particular observer.

    The *apparent position* of a body is its position relative to an
    observer adjusted for light-time delay, deflection (light rays
    bending as they pass large masses like the Sun or Jupiter), and
    aberration (light slanting because of the observer's motion through
    space).

    Included in aberration is the relativistic transformation that takes
    the position out of the BCRS centered on the solar system barycenter
    and into the reference frame of the observer.  In the case of an
    Earth observer, the transform takes the coordinate into the GCRS.

    """
    def altaz(self, temperature_C=None, pressure_mbar='standard'):
        """Compute (alt, az, distance) relative to the observer's horizon

        The altitude returned is an `Angle` in degrees above the
        horizon, while the azimuth is the compass direction in degrees
        with north being 0 degrees and east being 90 degrees.

        """
        try:
            if self.observer_data:
                elevation_m = self.observer_data.elevation_m
                R = self.observer_data.altaz_rotation
            else:
                # The old-fashioned way
                topos = self.observer.topos
                elevation_m = topos.elevation.m
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
                pressure_mbar = 1010.0 * exp(-elevation_m / 9.1e3)
            alt = refract(alt * RAD2DEG, temperature_C, pressure_mbar)
            alt = Angle(degrees=alt)

        return alt, Angle(radians=az), Distance(r_au)


class Geocentric(ICRF):
    """An (x,y,z) position measured from the geocenter."""


def ITRF_to_GCRS(t, rITRF):  # todo: velocity

    # Todo: wobble

    spin = rot_z(t.gast * tau / 24.0)
    position = einsum('ij...,j...->i...', spin, array(rITRF))
    return einsum('ij...,j...->i...', t.MT, position)
