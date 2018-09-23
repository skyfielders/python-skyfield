"""Classes representing different kinds of astronomical position."""

from numpy import arccos, array, clip, einsum, empty, exp, full, nan

from .constants import ANGVEL, DAY_S, DEG2RAD, RAD2DEG, tau
from .data.spice import inertial_frames
from .earthlib import compute_limb_angle, refract, reverse_terra
from .functions import dots, from_polar, length_of, rot_x, rot_z, to_polar
from .relativity import add_aberration, add_deflection
from .timelib import Time
from .units import Angle, Distance, Velocity, _interpret_angle

_ECLIPJ2000 = inertial_frames['ECLIPJ2000']
_GALACTIC = inertial_frames['GALACTIC']

def build_position(position_au, velocity_au_per_d=None, t=None,
                   center=None, target=None, observer_data=None):
    if center == 0:
        cls = Barycentric
    elif center == 399:
        cls = Geocentric
    elif observer_data is not None:
        cls = Geometric
    else:
        cls = ICRF
    return cls(position_au, velocity_au_per_d, t, center, target, observer_data)


class ICRF(object):
    """An (x, y, z) position and velocity oriented to the ICRF axes.

    The ICRF is a permanent coordinate system that has superseded the
    old series of equinox-based systems like B1900 and B1950.  Its axes
    are aligned with the axes of J2000 to within 0.02 arcseconds, which
    is tighter than the accuracy of J2000 itself.

    """
    def __init__(self, position_au, velocity_au_per_d=None, t=None,
                 center=None, target=None, observer_data=None):
        self.t = t
        self.position = Distance(position_au)
        if velocity_au_per_d is None:
            velocity_au_per_d = full(self.position.au.shape, nan)
        self.velocity = Velocity(velocity_au_per_d)
        # TODO: are center and target useful? Then why are they not
        # propagated down to Astrometric and Apparent positions?
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
        # TODO: set center and target of result
        p = self.position.au - body.position.au
        if self.velocity is None or body.velocity is None:
            v = None
        else:
            v = body.velocity.au_per_d - self.velocity.au_per_d
        return ICRF(p, v, self.t)

    def __getitem__(self, i):
        # TODO: grab other data too
        return type(self)(self.position.au[:,i], t=self.t[i])

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
        array([ 90.,   0.,  90., 180.])

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

    def cirs_xyz(self, epoch):
        """Compute cartesian CIRS coordinates at a given epoch (x, y, z).

        Calculate coordinates in the Celestial Intermediate Reference System
        (CIRS), a dynamical coordinate system referenced to the Celestial
        Intermediate Origin (CIO). As this is a dynamical system it must be
        calculated at a specific epoch.
        """
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

        vector = einsum('ij...,j...->i...', epoch.C, self.position.au)
        return Distance(vector)

    def cirs_radec(self, epoch):
        """Get spherical CIRS coordinates at a given epoch (ra, dec, distance).

        Calculate coordinates in the Celestial Intermediate Reference System
        (CIRS), a dynamical coordinate system referenced to the Celestial
        Intermediate Origin (CIO). As this is a dynamical system it must be
        calculated at a specific epoch.
        """
        r_au, dec, ra = to_polar(self.cirs_xyz(epoch).au)

        return (Angle(radians=ra, preference='hours'),
                Angle(radians=dec, signed=True),
                Distance(r_au))

    def ecliptic_xyz(self, epoch=None):
        """Compute J2000 ecliptic position vector (x, y, z).

        If you instead want the coordinates referenced to the dynamical
        system defined by the Earth's true equator and equinox, provide
        an epoch time.

        """
        if epoch is None:
            vector = _ECLIPJ2000.dot(self.position.au)
            return Distance(vector)

        position_au = self.position.au

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

        oblm, oblt, eqeq, psi, eps = epoch._earth_tilt
        e = oblt * DEG2RAD
        rotation = einsum('ij...,jk...->ik...', rot_x(-e), epoch.M)
        position_au = einsum('ij...,j...->i...', rotation, position_au)
        return Distance(position_au)

    def ecliptic_velocity(self):
        """Compute J2000 ecliptic velocity vector (x_dot, y_dot, z_dot)"""
        vector = _ECLIPJ2000.dot(self.velocity.au_per_d)
        return Velocity(vector)

    def ecliptic_latlon(self, epoch=None):
        """Compute J2000 ecliptic coordinates (lat, lon, distance)

        If you instead want the coordinates referenced to the dynamical
        system defined by the Earth's true equator and equinox, provide
        an epoch time.

        """
        vector = self.ecliptic_xyz(epoch)
        d, lat, lon = to_polar(vector.au)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def galactic_xyz(self):
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

    # Aliases; maybe someday turn into deprecations with warnings?
    ecliptic_position = ecliptic_xyz
    galactic_position = galactic_xyz

    def to_skycoord(self, unit=None):
        """Convert this distance to an AstroPy ``SkyCoord`` object."""
        from astropy.coordinates import SkyCoord
        from astropy.units import au
        x, y, z = self.position.au
        return SkyCoord(representation='cartesian', x=x, y=y, z=z, unit=au)

    def _to_spice_frame(self, name):
        vector = self.position.au
        vector = inertial_frames[name].dot(vector)
        d, dec, ra = to_polar(vector)
        return (Angle(radians=ra, preference='hours', signed=True),
                Angle(radians=dec),
                Distance(au=d))

    def from_altaz(self, alt=None, az=None, alt_degrees=None, az_degrees=None,
                   distance=Distance(au=0.1)):
        """Generate an Apparent position from an altitude and azimuth.

        The altitude and azimuth can each be provided as an `Angle`
        object, or else as a number of degrees provided as either a
        float or a tuple of degrees, arcminutes, and arcseconds::

            alt=Angle(...), az=Angle(...)
            alt_degrees=23.2289, az_degrees=142.1161
            alt_degrees=(23, 13, 44.1), az_degrees=(142, 6, 58.1)

        The distance should be a :class:`~skyfield.units.Distance`
        object, if provided; otherwise a default of 0.1 au is used.

        """
        # TODO: should this method live on another class?
        R = self.observer_data.altaz_rotation if self.observer_data else None
        if R is None:
            raise ValueError('only a position generated by a topos() call'
                             ' knows the orientation of the horizon'
                             ' and can understand altitude and azimuth')
        alt = _interpret_angle('alt', alt, alt_degrees)
        az = _interpret_angle('az', az, az_degrees)
        r = distance.au
        p = from_polar(r, alt, az)
        p = einsum('ji...,j...->i...', R, p)
        return Apparent(p)


# For compatibility with my original name for the class.  Not an
# important enough change to warrant a deprecation error for users, so:
ICRS = ICRF


class Geometric(ICRF):
    """An (x,y,z) vector between two instantaneous position.

    A geometric position is the difference between the Solar System
    positions of two bodies at exactly the same instant.  It is *not*
    corrected for the fact that, in real physics, it will take time for
    light to travel from one position to the other.

    """
    def altaz(self, temperature_C=None, pressure_mbar='standard'):
        """Compute (alt, az, distance) relative to the observer's horizon

        The altitude returned is an `Angle` in degrees above the
        horizon, while the azimuth is the compass direction in degrees
        with north being 0 degrees and east being 90 degrees.

        """
        return _to_altaz(self.position.au, self.observer_data,
                         temperature_C, pressure_mbar)


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
    <Barycentric position and velocity at date t center=0 target=499>

    """
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
        p, v, light_time = body._observe_from_bcrs(self)
        t = self.t
        astrometric = Astrometric(p, v, t, observer_data=self.observer_data)
        astrometric.light_time = light_time
        return astrometric

# TODO: pre-create a Barycentric object representing the SSB, and make
# it possible for it to observe() a planet.

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
        <Apparent position and velocity at date t>

        These transforms convert the position from the BCRS reference
        frame of the Solar System barycenter and to the reference frame
        of the observer.  In the specific case of an Earth observer, the
        output reference frame is the GCRS.

        """
        t = self.t
        position_au = self.position.au.copy()
        observer_data = self.observer_data
        gcrs_position = observer_data.gcrs_position

        if gcrs_position is None:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                position_au, gcrs_position)
            include_earth_deflection = nadir_angle >= 0.8

        add_deflection(position_au, observer_data.bcrs_position,
                       observer_data.ephemeris, t, include_earth_deflection)

        add_aberration(position_au, observer_data.bcrs_velocity,
                       self.light_time)

        return Apparent(position_au, t=t, observer_data=observer_data)


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
        return _to_altaz(self.position.au, self.observer_data,
                         temperature_C, pressure_mbar)


class Geocentric(ICRF):
    """An (x,y,z) position measured from the geocenter."""

    def subpoint(self):
        """Return the latitude and longitude directly beneath this position.

        Returns a :class:`~skyfield.toposlib.Topos` whose ``longitude``
        and ``latitude`` are those of the point on the Earth's surface
        directly beneath this position, and whose ``elevation`` is the
        height of this position above the Earth's surface.

        """
        if self.center != 399:  # TODO: should an __init__() check this?
            raise ValueError("you can only ask for the geographic subpoint"
                             " of a position measured from Earth's center")
        t = self.t
        xyz_au = einsum('ij...,j...->i...', t.M, self.position.au)
        lat, lon, elevation_m = reverse_terra(xyz_au, t.gast)

        # TODO. Move VectorFunction and Topos into this file, since the
        # three kinds of class work together: Topos is-a VF; VF.at() can
        # return a Geocentric position; and Geocentric.subpoint() should
        # return a Topos. I'm deferring the refactoring for now, to get
        # this new feature to users more quickly.
        from .toposlib import Topos
        return Topos(latitude=Angle(radians=lat),
                     longitude=Angle(radians=lon),
                     elevation_m=elevation_m)


def _to_altaz(position_au, observer_data, temperature_C, pressure_mbar):
    """Compute (alt, az, distance) relative to the observer's horizon.

    """
    elevation_m = observer_data.elevation_m
    R = observer_data.altaz_rotation

    if (elevation_m is None) or (R is None):
        raise ValueError('to compute an altazimuth position, you must'
                         ' observe from a specific Earth location that'
                         ' you specify using a Topos instance')

    # TODO: wobble

    position_au = einsum('ij...,j...->i...', R, position_au)
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

def ITRF_to_GCRS(t, rITRF):

    # Todo: wobble

    spin = rot_z(t.gast * tau / 24.0)
    position = einsum('ij...,j...->i...', spin, array(rITRF))
    return einsum('ij...,j...->i...', t.MT, position)

def ITRF_to_GCRS2(t, rITRF, vITRF):

    # Todo: wobble

    spin = rot_z(t.gast * tau / 24.0)

    position = einsum('ij...,j...->i...', spin, array(rITRF))
    position = einsum('ij...,j...->i...', t.MT, position)

    velocity = einsum('ij...,j...->i...', spin, array(vITRF))
    velocity = einsum('ij...,j...->i...', t.MT, velocity)
    velocity[0] += DAY_S * ANGVEL * - position[1]
    velocity[1] += DAY_S * ANGVEL * position[0]

    return position, velocity
