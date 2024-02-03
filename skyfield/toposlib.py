# -*- coding: utf-8 -*-
"""Classes that represent a ‘topocentric’ position on the Earth’s surface."""

from numpy import arctan2, array, array2string, cos, exp, sin, sqrt
from .constants import ANGVEL, DAY_S, RAD2DEG, T0, pi, tau
from .earthlib import refract
from .framelib import itrs
from .functions import (
    _T, angular_velocity_matrix, mxm, mxv, rot_y, rot_z,
)
from .descriptorlib import reify
from .units import Angle, Distance, _ltude
from .vectorlib import VectorFunction

_EARTH_ANGULAR_VELOCITY_VECTOR = array((0, 0, DAY_S * ANGVEL))

_lat_options = {'precision': 4, 'floatmode': 'fixed', 'sign': '+',
                'threshold': 5, 'edgeitems': 2}
_lon_options = {'precision': 4, 'floatmode': 'fixed',
                'threshold': 5, 'edgeitems': 2}
_elev_options = {'precision': 1, 'floatmode': 'fixed',
                 'threshold': 5, 'edgeitems': 2}

class ITRSPosition(VectorFunction):
    """An |xyz| position in the Earth-centered Earth-fixed (ECEF) ITRS frame."""

    center = 399

    def __init__(self, itrs_xyz):
        self.itrs_xyz = itrs_xyz
        x, y, z = itrs_xyz.au
        self._velocity_au_per_d = ANGVEL * DAY_S * array((-y, x, 0.0 * z))

    @property
    def target(self):
        # When used as a vector function, this Earth geographic location
        # computes positions from the Earth's center to itself.  (This
        # is a property, rather than an attribute, to avoid a circular
        # reference that delays garbage collection.)
        return self

    def _at(self, t):
        """Compute GCRS position and velocity at time `t`."""
        r = self.itrs_xyz.au
        v = self._velocity_au_per_d

        RT = _T(itrs.rotation_at(t))
        r = mxv(RT, r)
        v = mxv(RT, v)
        return r, v, None, None

class GeographicPosition(ITRSPosition):
    """A latitude-longitude-elevation position on Earth.

    Each instance of this class holds an |xyz| vector for a geographic
    position on, above, or below the Earth’s surface, in the ITRS
    reference frame: the international standard for an Earth-centered
    Earth-fixed (ECEF) reference frame.  Instead of instantiating this
    class directly, Skyfield users usually give a reference geoid the
    longitude and latitude they are interested in::

        from skyfield.api import wgs84
        topos = wgs84.latlon(37.3414, -121.6429)

    Once a geographic position has been created, here are its attributes
    and methods:

    """
    vector_name = 'Geodetic'

    def __init__(self, model, latitude, longitude, elevation, itrs_xyz):
        super(GeographicPosition, self).__init__(itrs_xyz)
        self.model = model
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self._R_lat = _R_lat = rot_y(latitude.radians)[::-1]
        self._R_latlon = mxm(_R_lat, rot_z(-longitude.radians))

    @reify  # not @property, so users have the option of overwriting it
    def target_name(self):
        return '{0} latitude {1} N longitude {2} E elevation {3} m'.format(
            self.model.name,
            array2string(self.latitude.degrees, **_lat_options),
            array2string(self.longitude.degrees, **_lon_options),
            array2string(self.elevation.m, **_elev_options))

    def lst_hours_at(self, t):
        """Return the Local Apparent Sidereal Time, in hours, at time ``t``.

        This location’s Local Apparent Sidereal Time (LAST) is the right
        ascension of the zenith at the time ``t``, as measured against
        the “true” Earth equator and equinox (rather than the fictional
        “mean” equator and equinox, which ignore the Earth’s nutation).

        """
        sprime = -47.0e-6 * (t.whole - T0 + t.tdb_fraction) / 36525.0
        return (t.gast + self.longitude._hours + sprime / 54000.0) % 24.0

    def refract(self, altitude_degrees, temperature_C, pressure_mbar):
        """Predict how the atmosphere will refract a position.

        Given a body that is standing ``altitude_degrees`` above the
        true horizon, return an ``Angle`` predicting its apparent
        altitude given the supplied temperature and pressure, either of
        which can be the string ``'standard'`` to use 10°C and a
        pressure of 1010 mbar adjusted for the elevation of this
        geographic location.

        """
        if temperature_C == 'standard':
            temperature_C = 10.0
        if pressure_mbar == 'standard':
            pressure_mbar = 1010.0 * exp(-self.elevation.m / 9.1e3)
        alt = refract(altitude_degrees, temperature_C, pressure_mbar)
        return Angle(degrees=alt)

    def rotation_at(self, t):
        """Compute rotation from GCRS to this location’s altazimuth system."""
        return mxm(self._R_latlon, itrs.rotation_at(t))

    def _dRdt_times_RT_at(self, t):
        # TODO: taking the derivative of the instantaneous angular
        # velocity would provide a more accurate transform.
        R = mxv(self._R_lat, _EARTH_ANGULAR_VELOCITY_VECTOR)
        return angular_velocity_matrix(R)

class Geoid(object):
    """An Earth ellipsoid: maps latitudes and longitudes to |xyz| positions.

    Instead of creating their own geoid object, most Skyfield users
    simply use the `wgs84` object that comes built-in.

    The math for turning a position into latitude and longitude is based
    on Dr. T.S. Kelso's quite helpful article `Orbital Coordinate
    Systems, Part III <https://www.celestrak.org/columns/v02n03/>`_.

    """
    def __init__(self, name, radius_m, inverse_flattening):
        self.name = name
        self.radius = Distance(m=radius_m)
        self.inverse_flattening = inverse_flattening
        omf = (inverse_flattening - 1.0) / inverse_flattening
        self._one_minus_flattening_squared = omf * omf
        f = 1.0 / inverse_flattening
        self._e2 = 2.0*f - f*f

    @reify
    def polar_radius(self):
        """The Earth’s polar radius, as a :class:`~skyfield.units.Distance`."""
        return Distance(self.radius.au * (1.0 - 1.0 / self.inverse_flattening))

    def latlon(self, latitude_degrees, longitude_degrees, elevation_m=0.0,
               cls=GeographicPosition):
        """Return a `GeographicPosition` for a given latitude and longitude.

        The longitude and latitude should both be specified in degrees.
        If no elevation in meters is supplied, the returned position
        will lie on the surface of the ellipsoid.  Longitude is positive
        towards the east, so supply a negative number for west::

            from skyfield.api import wgs84
            observatory = wgs84.latlon(37.3414, -121.6429)  # 121.6° West

        You can avoid remembering which directions are negative by using
        Skyfield’s compass direction constants, which have the values +1
        and −1::

            from skyfield.api import N, S, E, W
            observatory = wgs84.latlon(37.3414 * N, 121.6429 * W)

        """
        latitude = Angle(degrees=latitude_degrees)
        longitude = Angle(degrees=longitude_degrees)
        elevation = Distance(m=elevation_m)

        lat = latitude.radians
        lon = longitude.radians
        radius_au = self.radius.au
        elevation_au = elevation.au

        sinphi = sin(lat)
        cosphi = cos(lat)
        omf2 = self._one_minus_flattening_squared
        c = 1.0 / sqrt(cosphi * cosphi + sinphi * sinphi * omf2)
        s = omf2 * c

        # At equator: 6378 km, the Earth's actual radius at the equator.
        # At the pole: 6399 km, the Earth's radius of curvature at the pole.
        radius_xy = radius_au * c
        xy = (radius_xy + elevation_au) * cosphi
        x = xy * cos(lon)
        y = xy * sin(lon)

        # At equator: 6335 km, the Earth's radius of curvature at the equator.
        # At the pole: 6357 km, the Earth's actual radius at the pole.
        radius_z = radius_au * s
        z = (radius_z + elevation_au) * sinphi

        r = array((x, y, z))
        return cls(self, latitude, longitude, elevation, Distance(r))

    def latlon_of(self, position):
        """Return the latitude and longitude of a ``position``.

        The position’s ``.center`` must be 399, the center of the Earth.
        Geodetic latitude and longitude are returned as a pair of
        :class:`~skyfield.units.Angle` objects.

        """
        xyz_au, x, y, R, aC, hyp, lat = self._compute_latitude(position)
        lon = (arctan2(y, x) - pi) % tau - pi
        return Angle(radians=lat), Angle(radians=lon)

    def height_of(self, position):
        """Return the height above the Earth’s ellipsoid of a ``position``.

        The position’s ``.center`` must be 399, the center of the Earth.
        A :class:`~skyfield.units.Distance` is returned giving the
        position’s geodetic height above the Earth’s surface.

        """
        xyz_au, x, y, R, aC, hyp, lat = self._compute_latitude(position)
        height_au = sqrt(hyp * hyp + R * R) - aC
        return Distance(height_au)

    def geographic_position_of(self, position):
        """Return the `GeographicPosition` of a ``position``.

        The position’s ``.center`` must be 399, the center of the Earth.
        A `GeographicPosition` is returned giving the position’s
        geodetic ``latitude`` and ``longitude``, and an ``elevation``
        above or below the surface of the ellipsoid.

        """
        xyz_au, x, y, R, aC, hyp, lat = self._compute_latitude(position)
        lon = (arctan2(y, x) - pi) % tau - pi
        height_au = sqrt(hyp * hyp + R * R) - aC
        return GeographicPosition(
            latitude=Angle(radians=lat),
            longitude=Angle(radians=lon),
            elevation=Distance(height_au),
            itrs_xyz=Distance(xyz_au),
            model=self,
        )

    def subpoint_of(self, position):
        """Return the point on the ellipsoid directly below a ``position``.

        The position’s ``.center`` must be 399, the center of the Earth.
        Returns a `GeographicPosition` giving the geodetic ``latitude``
        and ``longitude`` that lie directly below the input position,
        and an ``elevation`` above the ellipsoid of zero.

        """
        xyz_au, x, y, R, aC, hyp, lat = self._compute_latitude(position)
        lon = (arctan2(y, x) - pi) % tau - pi
        return self.latlon(lat * RAD2DEG, lon * RAD2DEG)

    def _compute_latitude(self, position):
        if position.center != 399:
            raise ValueError(
                'you can only calculate a geographic position from a'
                ' position which is geocentric (center=399), but this'
                ' position has a center of {0}'.format(position.center)
            )
        xyz_au = position.frame_xyz(itrs).au
        x, y, z = xyz_au
        a = self.radius.au
        e2 = self._e2
        R = sqrt(x*x + y*y)
        lat = arctan2(z, R)
        for iteration in 0,1,2:
            sin_lat = sin(lat)
            e2_sin_lat = e2 * sin_lat
            # At 0°, aC = 6378 km, Earth's actual radius at the equator.
            # At 90°, aC = 6399 km, Earth's radius of curvature at the pole.
            aC = a / sqrt(1.0 - e2_sin_lat * sin_lat)
            hyp = z + aC * e2_sin_lat
            lat = arctan2(hyp, R)
        return xyz_au, x, y, R, aC, hyp, lat

    subpoint = geographic_position_of  # deprecated method name

wgs84 = Geoid('WGS84', 6378137.0, 298.257223563)
iers2010 = Geoid('IERS2010', 6378136.6, 298.25642)

wgs84.__doc__ = """World Geodetic System 1984 `Geoid`.

This is the standard geoid used by the GPS system,
and is likely the standard that’s intended
if you are supplied a latitude and longitude
that don’t specify an alternative geoid.

"""
iers2010.__doc__ = 'International Earth Rotation Service 2010 `Geoid`.'

# Compatibility with old versions of Skyfield:

class Topos(GeographicPosition):
    """Deprecated: use ``wgs84.latlon()`` or ``iers2010.latlon()`` instead."""

    def __init__(self, latitude=None, longitude=None, latitude_degrees=None,
                 longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):

        if latitude_degrees is not None:
            pass
        elif isinstance(latitude, Angle):
            latitude_degrees = latitude.degrees
        elif isinstance(latitude, (str, float, tuple)):
            latitude_degrees = _ltude(latitude, 'latitude', 'N', 'S')
        else:
            raise TypeError('please provide either latitude_degrees=<float>'
                            ' or latitude=<skyfield.units.Angle object>'
                            ' with north being positive')

        if longitude_degrees is not None:
            pass
        elif isinstance(longitude, Angle):
            longitude_degrees = longitude.degrees
        elif isinstance(longitude, (str, float, tuple)):
            longitude_degrees = _ltude(longitude, 'longitude', 'E', 'W')
        else:
            raise TypeError('please provide either longitude_degrees=<float>'
                            ' or longitude=<skyfield.units.Angle object>'
                            ' with east being positive')

        # Sneaky: the model thinks it's creating an object when really
        # it's just calling our superclass __init__() for us.  Alas, the
        # crimes committed to avoid duplicating code!  (This is actually
        # quite clean compared to the other alternatives I tried.)
        iers2010.latlon(latitude_degrees, longitude_degrees, elevation_m,
                        super(Topos, self).__init__)

        self.R_lat = self._R_lat  # On this old class, it was public.

    def itrf_xyz(self):
        return self.itrs_xyz
