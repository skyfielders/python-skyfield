# -*- coding: utf-8 -*-
"""Classes representing different kinds of astronomical position."""

from numpy import array, einsum, full, reshape, nan, nan_to_num
from . import framelib
from .constants import ANGVEL, AU_M, C, ERAD, DAY_S, RAD2DEG, tau
from .data.spice import inertial_frames
from .descriptorlib import reify
from .earthlib import compute_limb_angle
from .functions import (
    _T, _to_array, angle_between, from_spherical,
    length_of, mxm, mxv, rot_z, to_spherical,
)
from .geometry import intersect_line_and_sphere
from .relativity import add_aberration, add_deflection
from .timelib import Time
from .units import Angle, Distance, Velocity, _interpret_angle

_ECLIPJ2000 = inertial_frames['ECLIPJ2000']
_GIGAPARSEC_AU = 206264806247096.38  # 1e9 * 360 * 3600 / tau

def build_position(position_au, velocity_au_per_d=None, t=None,
                   center=None, target=None):
    if center == 0:
        cls = Barycentric
    elif center == 399:
        cls = Geocentric
    elif hasattr(center, 'rotation_at'):  # and thus deserves an altaz() method
        cls = Geometric
    else:
        cls = ICRF
    return cls(position_au, velocity_au_per_d, t, center, target)

def position_of_radec(ra_hours, dec_degrees, distance_au=_GIGAPARSEC_AU,
                      epoch=None, t=None, center=None, target=None):
    """Build a position object from a right ascension and declination.

    If a specific ``distance_au`` is not provided, Skyfield returns a
    position vector a gigaparsec in length.  This puts the position at a
    great enough distance that it will stand at the same right ascension
    and declination from any viewing position in the Solar System, to
    very high precision (within a few hundredths of a microarcsecond).

    If an ``epoch`` is specified, the input coordinates are understood
    to be in the dynamical system of that particular date.  Otherwise,
    they will be assumed to be ICRS (the modern replacement for J2000).

    """
    theta = _to_array(dec_degrees) / 360.0 * tau
    phi = _to_array(ra_hours) / 24.0 * tau
    position_au = from_spherical(distance_au, theta, phi)
    if epoch is not None:
        position_au = mxv(epoch.MT, position_au)
    return build_position(position_au, None, t, center, target)

def position_from_radec(ra_hours, dec_degrees, distance=1.0, epoch=None,
                        t=None, center=None, target=None):
    """DEPRECATED version of ``position_of_radec()``.

    Problems:

    * The ``distance`` parameter specifies no unit, contrary to Skyfield
      best practices.  I have no idea what I was thinking.

    * The default ``distance`` is far too small, since most objects for
      which users specify an RA and declination are out on the celestial
      sphere.  The hope was that users would see the length 1.0 and
      think, “ah, yes, that’s obviously a fake placeholder value.”  But
      it’s more likely that users will not even check the distance, or
      maybe not even realize that a distance is involved.

    """
    return position_of_radec(ra_hours, dec_degrees, distance, epoch,
                             t, center, target)

class ICRF(object):
    """An (x,y,z) position and velocity oriented to the ICRF axes.

    The International Coordinate Reference Frame (ICRF) is a permanent
    reference frame that is the replacement for J2000.  Their axes agree
    to within 0.02 arcseconds.  It also supersedes older equinox-based
    systems like B1900 and B1950.

    Each instance of this class provides a ``.position`` vector and a
    ``.velocity`` vector that specify (x,y,z) coordinates along the axes
    of the ICRF.  A specific time ``.t`` might be specified or might be
    ``None``.

    """
    center_barycentric = None
    _observer_gcrs_au = None
    _default_center = None
    _ephemeris = None  # cached so we can compute how light is deflected

    def __init__(self, position_au, velocity_au_per_d=None, t=None,
                 center=None, target=None):
        self.t = t
        self.position = Distance(position_au)
        if velocity_au_per_d is None:
            velocity_au_per_d = full(self.position.au.shape, nan)
        self.velocity = Velocity(velocity_au_per_d)
        self.center = self._default_center if center is None else center
        self.target = target
        if center == 0:
            self.center_barycentric = self

    @classmethod
    def from_radec(cls, ra_hours, dec_degrees,
                   distance_au=_GIGAPARSEC_AU, epoch=None):
        theta = _to_array(dec_degrees) / 360.0 * tau
        phi = _to_array(ra_hours) / 24.0 * tau
        position_au = from_spherical(distance_au, theta, phi)
        if epoch is not None:
            position_au = mxv(epoch.MT, position_au)
        return cls(position_au)

    @classmethod
    def from_time_and_frame_vectors(cls, t, frame, distance, velocity):
        """Constructor: build a position from two vectors in a reference frame.

        * ``t`` — The :class:`~skyfield.timelib.Time` of the position.
        * ``frame`` — A reference frame listed at `reference_frames`.
        * ``distance`` — A `Distance` x,y,z vector in the given frame.
        * ``velocity`` — A `Velocity` ẋ,ẏ,ż vector in the given frame.

        """
        r, v = distance.au, velocity.au_per_d
        at = getattr(frame, '_dRdt_times_RT_at', None)
        if at is not None:
            V = at(t)
            v = v - mxv(V, r)  # subtract instead of transposing
        RT = _T(frame.rotation_at(t))
        r = mxv(RT, r)
        v = mxv(RT, v)
        return cls(r, v, t)  # TODO: args for center and target?

    def __repr__(self):
        name = self.__class__.__name__
        center = self.center
        if name == 'Barycentric' and center == 0:
            suffix = ' BCRS'
        elif name == 'Apparent' and center == 399:
            suffix = ' GCRS'
        elif name != 'ICRF':
            suffix = ' ICRS'
        else:
            suffix = ''

        center = self.center
        target = self.target

        center_name = getattr(center, 'center_name', None)
        if center_name is None:
            center_name = str(center)

        target_name = getattr(target, 'target_name', None)
        if target_name is None:
            target_name = str(target)

        return '<{0}{1} position{2}{3}{4}{5}>'.format(
            name,
            suffix,
            '' if (self.velocity is None) else ' and velocity',
            '' if self.t is None else ' at date t',
            '' if self.center is None else ' center={0}'.format(center_name),
            '' if self.target is None else ' target={0}'.format(target_name),
        )

    def __sub__(self, body):
        """Subtract two ICRF vectors to produce a third."""
        # TODO: set center and target of result
        p = self.position.au - body.position.au
        if self.velocity is None or body.velocity is None:
            v = None
        else:
            v = self.velocity.au_per_d - body.velocity.au_per_d
        return ICRF(p, v, self.t)

    def __getitem__(self, i):
        return type(self)(
            self.position.au[:,i],
            self.velocity.au_per_d[:,i],
            self.t[i],
            self.center,
            self.target,
        )

    def __neg__(self):
        return type(self)(
            -self.position.au,
            -self.velocity.au_per_d,
            self.t,
            self.target,
            self.center,
        )

    def distance(self):
        """Compute the distance from the origin to this position.

        The return value is a :class:`~skyfield.units.Distance` that
        prints itself out in astronomical units (au) but that also
        offers attributes ``au``, ``km``, and ``m`` if you want to
        access its magnitude as a number.

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

    @reify
    def light_time(self):
        """Length of this vector in days of light travel time."""
        return self.distance().m / C * DAY_S

    def radec(self, epoch=None):
        r"""Compute equatorial (RA, declination, distance)

        When called without a parameter, this returns standard ICRF
        right ascension and declination:

        >>> from skyfield.api import load
        >>> ts = load.timescale()
        >>> t = ts.utc(2020, 5, 13, 10, 32)
        >>> eph = load('de421.bsp')
        >>> astrometric = eph['earth'].at(t).observe(eph['sun'])
        >>> ra, dec, distance = astrometric.radec()
        >>> print(ra, dec, sep='\n')
        03h 21m 47.67s
        +18deg 28' 55.3"

        If you instead want the coordinates referenced to the dynamical
        system defined by the Earth's true equator and equinox, provide
        a specific epoch time.

        >>> ra, dec, distance = astrometric.apparent().radec(epoch='date')
        >>> print(ra, dec, sep='\n')
        03h 22m 54.73s
        +18deg 33' 04.5"

        To get J2000.0 coordinates, simply pass ``ts.J2000``.

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
            position_au = mxv(epoch.M, position_au)
        r_au, dec, ra = to_spherical(position_au)
        return (Angle(radians=ra, preference='hours'),
                Angle(radians=dec, signed=True),
                Distance(r_au))

    def separation_from(self, another_icrf):
        """Return the angle between this position and another.

        >>> from skyfield.api import load
        >>> ts = load.timescale()
        >>> t = ts.utc(2020, 4, 18)
        >>> eph = load('de421.bsp')
        >>> sun, venus, earth = eph['sun'], eph['venus'], eph['earth']
        >>> e = earth.at(t)
        >>> s = e.observe(sun)
        >>> v = e.observe(venus)
        >>> print(s.separation_from(v))
        43deg 23' 23.1"

        You can also compute separations across an array of positions.

        >>> t = ts.utc(2020, 4, [18, 19, 20])
        >>> e = earth.at(t)
        >>> print(e.observe(sun).separation_from(e.observe(venus)))
        3 values from 43deg 23' 23.1" to 42deg 49' 46.6"

        """
        u = self.position.au
        v = another_icrf.position.au

        # Allow an array of positions to be compared with a single other
        # position.
        difference = len(u.shape) - len(v.shape)
        if difference:
            if difference > 0:
                v = reshape(v, v.shape + (1,) * difference)
            else:
                u = reshape(u, u.shape + (1,) * -difference)

        return Angle(radians=angle_between(u, v))

    # TODO: build a reference frame for the following two methods.

    def cirs_xyz(self, epoch):
        """Compute cartesian CIRS coordinates at a given epoch (x,y,z).

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

        vector = mxv(epoch.C, self.position.au)
        return Distance(vector)

    def cirs_radec(self, epoch):
        """Get spherical CIRS coordinates at a given epoch (ra, dec, distance).

        Calculate coordinates in the Celestial Intermediate Reference System
        (CIRS), a dynamical coordinate system referenced to the Celestial
        Intermediate Origin (CIO). As this is a dynamical system it must be
        calculated at a specific epoch.
        """
        r_au, dec, ra = to_spherical(self.cirs_xyz(epoch).au)

        return (Angle(radians=ra, preference='hours'),
                Angle(radians=dec, signed=True),
                Distance(r_au))

    # Deprecated methods, that have been replaced by `framelib.py` plus
    # the "frame" methods in the next section.

    def ecliptic_xyz(self, epoch=None):
        if epoch is None:
            return self.frame_xyz(framelib.ecliptic_J2000_frame)
        return _Fake(self, epoch).frame_xyz(framelib.ecliptic_frame)
    def ecliptic_velocity(self):
        return self.frame_xyz_and_velocity(framelib.ecliptic_J2000_frame)[1]
    def ecliptic_latlon(self, epoch=None):
        if epoch is None:
            return self.frame_latlon(framelib.ecliptic_J2000_frame)
        return _Fake(self, epoch).frame_latlon(framelib.ecliptic_frame)

    def galactic_xyz(self): return self.frame_xyz(framelib.galactic_frame)
    def galactic_latlon(self): return self.frame_latlon(framelib.galactic_frame)
    ecliptic_position = ecliptic_xyz  # old alias
    galactic_position = galactic_xyz  # old alias

    # New methods for converting to and from `framelib.py` reference frames.

    def frame_xyz(self, frame):
        """Return this position as an (x,y,z) vector in a reference frame.

        Returns a :class:`~skyfield.units.Distance` object giving the
        (x,y,z) of this position in the given ``frame``.  See
        `reference_frames`.

        """
        return Distance(mxv(frame.rotation_at(self.t), self.position.au))

    def frame_xyz_and_velocity(self, frame):
        """Return (x,y,z) position and velocity vectors in a reference frame.

        Returns two vectors in the given coordinate ``frame``: a
        :class:`~skyfield.units.Distance` providing an (x,y,z) position
        and a :class:`~skyfield.units.Velocity` giving (xdot,ydot,zdot)
        velocity.  See `reference_frames`.

        """
        R = frame.rotation_at(self.t)
        r, v = self.position.au, self.velocity.au_per_d
        r = mxv(R, r)
        v = mxv(R, v)
        at = getattr(frame, '_dRdt_times_RT_at', None)
        if at is not None:
            V = at(self.t)
            v += mxv(V, r)
        return Distance(r), Velocity(v)

    def frame_latlon(self, frame):
        """Return longitude, latitude, and distance in the given frame.

        Returns a 3-element tuple giving the latitude and longitude as
        :class:`~skyfield.units.Angle` objects and the range to the
        target as a :class:`~skyfield.units.Distance`.  See
        `reference_frames`.

        """
        vector = mxv(frame.rotation_at(self.t), self.position.au)
        d, lat, lon = to_spherical(vector)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(d))

    def to_skycoord(self, unit=None):
        """Convert this distance to an AstroPy ``SkyCoord`` object."""
        from astropy.coordinates import SkyCoord
        from astropy.units import au
        x, y, z = self.position.au
        return SkyCoord(representation_type='cartesian', x=x, y=y, z=z, unit=au)

    def _to_spice_frame(self, name):
        vector = self.position.au
        vector = inertial_frames[name].dot(vector)
        d, dec, ra = to_spherical(vector)
        return (Angle(radians=ra, preference='hours', signed=True),
                Angle(radians=dec),
                Distance(au=d))

    def is_sunlit(self, ephemeris):
        """Return whether a position in Earth orbit is in sunlight.

        Returns ``True`` or ``False``, or an array of such values, to
        indicate whether this position is in sunlight or is blocked by
        the Earth’s shadow.  It should work with positions produced
        either by calling ``at()`` on a satellite object, or by calling
        ``at()`` on the relative position ``sat - topos`` of a satellite
        with respect to an Earth observer’s position.  See
        :ref:`satellite-is-sunlit`.

        """
        if self.center == 399:
            earth_m = - self.position.m
        else:
            gcrs_position = self._observer_gcrs_au
            if gcrs_position is None:
                raise ValueError('cannot tell whether this position is sunlit')
            earth_m = - self.position.m - gcrs_position * AU_M

        sun_m = (ephemeris['sun'] - ephemeris['earth']).at(self.t).position.m
        near, far = intersect_line_and_sphere(sun_m + earth_m, earth_m, ERAD)
        return nan_to_num(far) <= 0

    def is_behind_earth(self):
        """Return whether the Earth blocks the view of this object.

        For a position centered on an Earth-orbiting satellite, return
        whether the target is in eclipse behind the disc of the Earth.
        See :ref:`is-behind-earth`.

        """
        observer_gcrs_au = self._observer_gcrs_au
        if observer_gcrs_au is None:
            raise ValueError('can only compute Earth occultation for'
                             ' positions observed from an Earth satellite')
        earth_m = - observer_gcrs_au * AU_M
        vector_m = self.position.m
        near, far = intersect_line_and_sphere(vector_m, earth_m, ERAD)
        return nan_to_num(far) > 0

    @reify
    def _altaz_rotation(self):
        # Return and cache (with @reify) the orientation of this
        # observer, in case a single observer.at() position is used in
        # several subsequent .observe().apparent().altaz() calls.
        rotation_at = getattr(self.target, 'rotation_at', None)
        if rotation_at is None:
            raise ValueError(_altaz_message)
        return rotation_at(self.t)

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

        rotation_at = getattr(self.target, 'rotation_at', None)
        if rotation_at is None:
            raise ValueError(_altaz_message)
        R = rotation_at(self.t)

        alt = _interpret_angle('alt', alt, alt_degrees)
        az = _interpret_angle('az', az, az_degrees)
        r = distance.au
        p = from_spherical(r, alt, az)
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

    Both the ``.position`` and ``.velocity`` are (x,y,z) vectors
    oriented along the axes of the International Celestial Reference
    System (ICRS), the modern replacement for J2000 coordinates.

    """
    def altaz(self, temperature_C=None, pressure_mbar='standard'):
        """Compute (alt, az, distance) relative to the observer's horizon

        The altitude returned is an :class:`~skyfield.units.Angle`
        measured in degrees above the horizon, while the azimuth
        :class:`~skyfield.units.Angle` measures east along the horizon
        from geographic north (so 0 degrees means north, 90 is east, 180
        is south, and 270 is west).

        By default, Skyfield does not adjust the altitude for
        atmospheric refraction.  If you want Skyfield to estimate how
        high the atmosphere might lift the body's image, give the
        argument ``temperature_C`` either the temperature in degrees
        centigrade, or the string ``'standard'`` (in which case 10°C is
        used).

        When calculating refraction, Skyfield uses the observer’s
        elevation above sea level to estimate the atmospheric pressure.
        If you want to override that value, simply provide a number
        through the ``pressure_mbar`` parameter.

        """
        return _to_altaz(self, temperature_C, pressure_mbar)


class Barycentric(ICRF):
    """An (x,y,z) position measured from the Solar System barycenter.

    Skyfield generates a `Barycentric` position measured from the
    gravitational center of the Solar System whenever you ask a body for
    its location at a particular time:

    >>> t = ts.utc(2003, 8, 29)
    >>> mars.at(t)
    <Barycentric BCRS position and velocity at date t center=0 target=499>

    This class’s ``.position`` and ``.velocity`` are (x,y,z) vectors in
    the Barycentric Celestial Reference System (BCRS), the modern
    replacement for J2000 coordinates measured from the Solar System
    Barycenter.

    """
    _default_center = 0

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
        <Astrometric ICRS position and velocity at date t center=399 target=499>

        """
        p, v, t, light_time = body._observe_from_bcrs(self)
        astrometric = Astrometric(p, v, t, self.target, body.target)
        astrometric._ephemeris = self._ephemeris
        astrometric.center_barycentric = self
        astrometric.light_time = light_time
        return astrometric

# TODO: pre-create a Barycentric object representing the SSB, and make
# it possible for it to observe() a planet.

class Astrometric(ICRF):
    """An astrometric (x,y,z) position relative to a particular observer.

    The astrometric position of a body is its position relative to an
    observer, adjusted for light-time delay.  It is the position of the
    body back when it emitted (or reflected) the light that is now
    reaching the observer's eye or telescope.  Astrometric positions are
    usually generated in Skyfield by calling the `Barycentric` method
    `observe()`, which performs the light-time correction.

    Both the ``.position`` and ``.velocity`` are ``[x y z]`` vectors
    oriented along the axes of the ICRF, the modern replacement for the
    J2000 reference frame.

    It is common to either call ``.radec()`` (with no argument) on an
    astrometric position to generate an *astrometric place* right
    ascension and declination with respect to the ICRF axes, or else to
    call ``.apparent()`` to generate an :class:`Apparent` position.

    """
    def apparent(self):
        """Compute an :class:`Apparent` position for this body.

        This applies two effects to the position that arise from
        relativity and shift slightly where the other body will appear
        in the sky: the deflection that the image will experience if its
        light passes close to large masses in the Solar System, and the
        aberration of light caused by the observer's own velocity.

        >>> earth.at(t).observe(mars).apparent()
        <Apparent GCRS position and velocity at date t center=399 target=499>

        These transforms convert the position from the BCRS reference
        frame of the Solar System barycenter and to the reference frame
        of the observer.  In the specific case of an Earth observer, the
        output reference frame is the GCRS.

        """
        t = self.t
        target_au = self.position.au.copy()

        cb = self.center_barycentric
        bcrs_position = cb.position.au
        bcrs_velocity = cb.velocity.au_per_d
        observer_gcrs_au = cb._observer_gcrs_au

        # If a single observer position (3,) is observing an array of
        # targets (3,n), then deflection and aberration will complain
        # that "operands could not be broadcast together" unless we give
        # the observer another dimension too.
        if len(bcrs_position.shape) < len(target_au.shape):
            shape = bcrs_position.shape + (1,)
            bcrs_position = bcrs_position.reshape(shape)
            bcrs_velocity = bcrs_velocity.reshape(shape)
            if observer_gcrs_au is not None:
                observer_gcrs_au = observer_gcrs_au.reshape(shape)

        if observer_gcrs_au is None:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                target_au, observer_gcrs_au)
            include_earth_deflection = nadir_angle >= 0.8

        add_deflection(target_au, bcrs_position,
                       self._ephemeris, t, include_earth_deflection)

        add_aberration(target_au, bcrs_velocity, self.light_time)

        apparent = Apparent(target_au, None, t, self.center, self.target)
        apparent.center_barycentric = self.center_barycentric
        apparent._observer_gcrs_au = observer_gcrs_au
        return apparent

class Apparent(ICRF):
    """An apparent ``[x y z]`` position relative to a particular observer.

    This class’s vectors provide the position and velocity of a body
    relative to an observer, adjusted to predict where the body’s image
    will really appear (hence "apparent") in the sky:

    * Light-time delay, as already present in an `Astrometric` position.

    * Deflection: gravity bends light, and thus the image of a distant
      object, as the light passes massive objects like Jupiter, Saturn,
      and the Sun.  For an observer on the Earth’s surface or in Earth
      orbit, the slight deflection by the gravity of the Earth itself is
      also included.

    * Aberration: incoming light arrives slanted because of the
      observer's motion through space.

    These positions are usually produced in Skyfield by calling the
    `apparent()` method of an `Astrometric` object.

    Both the ``.position`` and ``.velocity`` are ``[x y z]`` vectors
    oriented along the axes of the ICRF, the modern replacement for the
    J2000 reference frame.  If the observer is at the geocenter, they
    are more specifically GCRS coordinates.  Two common coordinates that
    this vector can generate are:

    * *Proper place:* call ``.radec()`` without arguments to compute
      right ascension and declination with respect to the fixed axes of
      the ICRF.

    * *Apparent place,* the most popular option: call ``.radec('date')``
      to generate right ascension and declination with respect to the
      equator and equinox of date.

    """
    def altaz(self, temperature_C=None, pressure_mbar='standard'):
        """Compute (alt, az, distance) relative to the observer's horizon

        The altitude returned is an :class:`~skyfield.units.Angle`
        measured in degrees above the horizon, while the azimuth
        :class:`~skyfield.units.Angle` measures east along the horizon
        from geographic north (so 0 degrees means north, 90 is east, 180
        is south, and 270 is west).

        By default, Skyfield does not adjust the altitude for
        atmospheric refraction.  If you want Skyfield to estimate how
        high the atmosphere might lift the body's image, give the
        argument ``temperature_C`` either the temperature in degrees
        centigrade, or the string ``'standard'`` (in which case 10°C is
        used).

        When calculating refraction, Skyfield uses the observer’s
        elevation above sea level to estimate the atmospheric pressure.
        If you want to override that value, simply provide a number
        through the ``pressure_mbar`` parameter.

        """
        return _to_altaz(self, temperature_C, pressure_mbar)


class Geocentric(ICRF):
    """An (x,y,z) position measured from the center of the Earth.

    A geocentric position is the difference between the position of the
    Earth at a given instant and the position of a target body at the
    same instant, without accounting for light-travel time or the effect
    of relativity on the light itself.

    Its ``.position`` and ``.velocity`` vectors have (x,y,z) axes that
    are those of the Geocentric Celestial Reference System (GCRS), an
    inertial system that is an update to J2000 and that does not rotate
    with the Earth itself.

    """
    _default_center = 399

    def itrf_xyz(self):
        """Deprecated; instead, call ``.frame_xyz(itrs)``. \
        See `reference_frames`."""
        return self.frame_xyz(framelib.itrs)

    def subpoint(self):
        """Deprecated; instead, call either ``iers2010.subpoint(pos)`` or \
        ``wgs84.subpoint(pos)``."""
        from .toposlib import iers2010
        return iers2010.subpoint(self)

def _to_altaz(position, temperature_C, pressure_mbar):
    """Compute (alt, az, distance) relative to the observer's horizon."""
    cb = position.center_barycentric
    if cb is not None:
        R = cb._altaz_rotation
    else:
        rotation_at = getattr(position.center, 'rotation_at')
        if rotation_at is not None:
            R = rotation_at(position.t)
        else:
            raise ValueError(_altaz_message)

    position_au = mxv(R, position.position.au)
    r_au, alt, az = to_spherical(position_au)

    if temperature_C is None:
        alt = Angle(radians=alt)
    else:
        refract = getattr(position.center, 'refract', None)
        if refract is None:
            raise ValueError(_altaz_message)
        alt = position.center.refract(
            alt * RAD2DEG, temperature_C, pressure_mbar,
        )

    return alt, Angle(radians=az), Distance(r_au)

_altaz_message = (
    'to compute an altazimuth position, you must observe from a'
    ' specific Earth location or from a position on another body'
    ' loaded from a set of planetary constants'
)

class _Fake(ICRF):  # support for deprecated frame rotation methods above
    def __init__(self, original, epoch):
        self.position = original.position
        if isinstance(epoch, Time):
            self.t = epoch
        elif isinstance(epoch, float):
            self.t = Time(None, tt=epoch)
        elif epoch == 'date':
            self.t = original.t
        else:
            raise ValueError('the epoch= must be a Time object,'
                             ' a floating point Terrestrial Time (TT),'
                             ' or the string "date" for epoch-of-date')

def ITRF_to_GCRS(t, rITRF):  # Deprecated; for compatibility with old versions.
    return mxv(_T(framelib.itrs.rotation_at(t)), rITRF)

def ITRF_to_GCRS2(t, rITRF, vITRF, _high_accuracy=False):
    position = array(rITRF)
    velocity = array(vITRF)

    # TODO: This is expensive, and should be extensively trimmed to only
    # include the most important terms underlying GAST.  But it improves
    # the precision by something like 1e5 times when compared to using
    # the round number skyfield.constants.ANGVEL!
    #
    # See the test `test_velocity_in_ITRF_to_GCRS2()`.
    #
    if _high_accuracy:
        _one_second = 1.0 / DAY_S
        t_later = t.ts.tt_jd(t.whole, t.tt_fraction + _one_second)
        angvel = (t_later.gast - t.gast) / 24.0 * tau
    else:
        angvel = ANGVEL

    spin = rot_z(t.gast / 24.0 * tau)
    R = mxm(t.MT, spin)

    z = 0.0 * angvel
    import numpy as np
    V = np.array((
        (z,-DAY_S * angvel,z),
        (DAY_S * angvel,z,z),
        (z,z,z),
    ))

    velocity = velocity + mxv(V, position)
    position = mxv(R, position)
    velocity = mxv(R, velocity)

    return position, velocity
