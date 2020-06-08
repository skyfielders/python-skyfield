# -*- coding: utf-8 -*-
"""Classes representing different kinds of astronomical position."""

from numpy import arccos, array, clip, einsum, exp, full, nan, nan_to_num
from .constants import ANGVEL, AU_M, ERAD, DAY_S, RAD2DEG, tau
from .data.spice import inertial_frames
from .earthlib import compute_limb_angle, refract, reverse_terra
from .geometry import intersect_line_and_sphere
from .functions import (
    _mxv, _mxm, dots, from_spherical, length_of, rot_x, rot_z,
    to_spherical, _to_array,
)
from .relativity import add_aberration, add_deflection
from .timelib import Time
from .units import Angle, Distance, Velocity, _interpret_angle

_ECLIPJ2000 = inertial_frames['ECLIPJ2000']
_GALACTIC = inertial_frames['GALACTIC']
_GIGAPARSEC_AU = 206264806247096.38  # 1e9 * 360 * 3600 / tau

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

def position_of_radec(ra_hours, dec_degrees, distance_au=_GIGAPARSEC_AU,
                      epoch=None, t=None, center=None, target=None,
                      observer_data=None):
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
        position_au = einsum('ij...,j...->i...', epoch.MT, position_au)
    return build_position(position_au, None, t, center, target, observer_data)

def position_from_radec(ra_hours, dec_degrees, distance=1.0, epoch=None,
                        t=None, center=None, target=None, observer_data=None):
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
                             t, center, target, observer_data)

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
    _default_center = None

    def __init__(self, position_au, velocity_au_per_d=None, t=None,
                 center=None, target=None, observer_data=None):
        self.t = t
        self.position = Distance(position_au)
        if velocity_au_per_d is None:
            velocity_au_per_d = full(self.position.au.shape, nan)
        self.velocity = Velocity(velocity_au_per_d)
        # TODO: are center and target useful? Then why are they not
        # propagated down to Astrometric and Apparent positions?
        self.center = self._default_center if center is None else center
        self.target = target
        self.observer_data = observer_data

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
        return '<{0}{1} position{2}{3}{4}{5}>'.format(
            name,
            suffix,
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
            v = self.velocity.au_per_d - body.velocity.au_per_d
        return ICRF(p, v, self.t)

    def __getitem__(self, i):
        return type(self)(
            self.position.au[:,i],
            self.velocity.au_per_d[:,i],
            self.t[i],
            self.center,
            self.target,
            # TODO: figure out how to slice observer data
            #self.observer_data,
        )

    def __neg__(self):
        return type(self)(
            -self.position.au,
            -self.velocity.au_per_d,
            self.t,
            self.target,
            self.center,
            # TODO: figure out how to invert observer data
            #self.observer_data,
        )

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

        >>> from skyfield.api import load
        >>> ts = load.timescale(builtin=True)
        >>> t = ts.utc(2020, 5, 13, 10, 32)
        >>> eph = load('de421.bsp')
        >>> astrometric = eph['earth'].at(t).observe(eph['sun'])
        >>> ra, dec, distance = astrometric.radec()
        >>> print(ra, dec, sep='\n')
        03h 21m 47.67s
        +18deg 28' 55.3"

        If you instead want the coordinates referenced to the dynamical
        system defined by the Earth's mean equator and equinox, provide
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
            position_au = einsum('ij...,j...->i...', epoch.M, position_au)
        r_au, dec, ra = to_spherical(position_au)
        return (Angle(radians=ra, preference='hours'),
                Angle(radians=dec, signed=True),
                Distance(r_au))

    def separation_from(self, another_icrf):
        """Return the angle between this position and another.

        >>> from skyfield.api import load
        >>> ts = load.timescale(builtin=True)
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

        vector = einsum('ij...,j...->i...', epoch.C, self.position.au)
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

    def ecliptic_xyz(self, epoch=None):
        """Compute J2000 ecliptic position vector (x,y,z).

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

        _, d_eps = epoch._nutation_angles_radians
        true_obliquity = epoch._mean_obliquity_radians + d_eps
        rotation = _mxm(rot_x(- true_obliquity), epoch.M)
        position_au = _mxv(rotation, position_au)
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
        d, lat, lon = to_spherical(vector.au)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def galactic_xyz(self):
        """Compute galactic coordinates (x,y,z)"""
        vector = _GALACTIC.dot(self.position.au)
        return Distance(vector)

    def galactic_latlon(self):
        """Compute galactic coordinates (lat, lon, distance)"""
        vector = _GALACTIC.dot(self.position.au)
        d, lat, lon = to_spherical(vector)
        return (Angle(radians=lat, signed=True),
                Angle(radians=lon),
                Distance(au=d))

    def frame_xyz(self, frame):
        """Express this position as an (x,y,z) vector in a particular frame."""
        R = frame.rotation_at(self.t)
        return Distance(au=_mxv(R, self.position.au))

    def frame_latlon(self, frame):
        """Return as longitude, latitude, and distance in the given frame."""
        R = frame.rotation_at(self.t)
        vector = _mxv(R, self.position.au)
        d, lat, lon = to_spherical(vector)
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
            o = self.observer_data
            if (o is not None) and (o.gcrs_position is not None):
                earth_m = - self.position.m - o.gcrs_position * AU_M
            else:
                raise ValueError('cannot tell whether this position is sunlit')

        sun_m = (ephemeris['sun'] - ephemeris['earth']).at(self.t).position.m
        near, far = intersect_line_and_sphere(sun_m + earth_m, earth_m, ERAD)
        return nan_to_num(far) <= 0

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
    oriented along the axes of the International Terrestrial Reference
    Frame (ITRF), the modern replacement for J2000 coordinates.

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
        return _to_altaz(self.position.au, self.observer_data,
                         temperature_C, pressure_mbar)


class Barycentric(ICRF):
    """An (x,y,z) position measured from the Solar System barycenter.

    Skyfield generates a `Barycentric` position measured from the
    gravitational center of the Solar System whenever you ask a body for
    its location at a particular time:

    >>> t = ts.utc(2003, 8, 29)
    >>> mars.at(t)
    <Barycentric BCRS position and velocity at date t center=0 target=499>

    Both the ``.position`` and ``.velocity`` are (x,y,z) vectors
    oriented along the axes of the International Terrestrial Reference
    Frame (ITRF), the modern replacement for J2000 coordinates.

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
        <Astrometric ICRS position and velocity at date t center=399 target=499>

        """
        p, v, t, light_time = body._observe_from_bcrs(self)
        astrometric = Astrometric(p, v, t,
                                  center=self.target, target=body.target,
                                  observer_data=self.observer_data)
        astrometric.light_time = light_time
        return astrometric

# TODO: pre-create a Barycentric object representing the SSB, and make
# it possible for it to observe() a planet.

class Astrometric(ICRF):
    """An astrometric (x,y,z) position relative to a particular observer.

    The *astrometric position* of a body is its position relative to an
    observer, adjusted for light-time delay.  It is the position of the
    body back when it emitted (or reflected) the light that is now
    reaching the observer's eyes or telescope.

    Both the ``.position`` and ``.velocity`` are (x,y,z) vectors
    oriented along the axes of the International Terrestrial Reference
    Frame (ITRF), the modern replacement for J2000 coordinates.

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
        <Apparent GCRS position and velocity at date t center=399 target=499>

        These transforms convert the position from the BCRS reference
        frame of the Solar System barycenter and to the reference frame
        of the observer.  In the specific case of an Earth observer, the
        output reference frame is the GCRS.

        """
        t = self.t
        target_au = self.position.au.copy()

        observer_data = self.observer_data
        gcrs_position = observer_data.gcrs_position
        bcrs_position = observer_data.bcrs_position
        bcrs_velocity = observer_data.bcrs_velocity

        # If a single observer position (3,) is observing an array of
        # targets (3,n), then deflection and aberration will complain
        # that "operands could not be broadcast together" unless we give
        # the observer another dimension too.
        if len(bcrs_position.shape) < len(target_au.shape):
            shape = bcrs_position.shape + (1,)
            bcrs_position = bcrs_position.reshape(shape)
            bcrs_velocity = bcrs_velocity.reshape(shape)
            if gcrs_position is not None:
                gcrs_position = gcrs_position.reshape(shape)

        if gcrs_position is None:
            include_earth_deflection = array((False,))
        else:
            limb_angle, nadir_angle = compute_limb_angle(
                target_au, gcrs_position)
            include_earth_deflection = nadir_angle >= 0.8

        add_deflection(target_au, bcrs_position,
                       observer_data.ephemeris, t, include_earth_deflection)

        add_aberration(target_au, bcrs_velocity, self.light_time)

        return Apparent(target_au, t=t,
                        center=self.center, target=self.target,
                        observer_data=observer_data)


class Apparent(ICRF):
    """An apparent (x,y,z) position relative to a particular observer.

    The *apparent position* of a body is its position relative to an
    observer adjusted for light-time delay, deflection (light rays
    bending as they pass large masses like the Sun or Jupiter), and
    aberration (light slanting because of the observer's motion through
    space).

    Included in aberration is the relativistic transformation that takes
    the ``.position`` and ``.velocity`` (x,y,z) vectors out of the BCRS,
    centered on the Solar System barycenter, and into the reference
    frame of the observer.  In the case of an Earth observer, the
    transform takes the coordinate into the GCRS.

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
        return _to_altaz(self.position.au, self.observer_data,
                         temperature_C, pressure_mbar)


class Geocentric(ICRF):
    """An (x,y,z) position measured from the center of the Earth.

    A geocentric position is the difference between the position of the
    Earth at a given instant and the position of a target body at the
    same instant, without accounting for light-travel time or the effect
    of relativity on the light itself.

    Its ``.position`` and ``.velocity`` vectors have (x,y,z) axes that
    are those of the International Terrestrial Reference Frame (ITRF),
    an inertial system that is an update to J2000 and that does not
    rotate with the Earth itself.

    """
    _default_center = 399

    def itrf_xyz(self):
        """Return this position as an (x,y,z) vector in the ITRF frame.

        Returns a :class:`~skyfield.units.Distance` object giving the
        (x,y,z) of this coordinate in the International Terrestrial
        Reference Frame (ITRF), an internationally agreed upon
        Earth-centered Earth-fixed (ECEF) coordinate system that
        rotates with the surface of the Earth itself.

        """
        if self.center != 399:
            raise ValueError("you can only ask for an Earth-centered position"
                             " to be converted into an ITRF coordinate")

        t = self.t
        au = einsum('ij...,j...->i...', t.M, self.position.au)

        spin = rot_z(- t.gast * tau / 24.0)
        au = einsum('ij...,j...->i...', spin, array(au))

        return Distance(au)

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
                         ' observe from a specific Earth location or from'
                         ' a position on another body loaded from a set'
                         ' of planetary constants')

    # TODO: wobble

    position_au = einsum('ij...,j...->i...', R, position_au)
    r_au, alt, az = to_spherical(position_au)

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
