# -*- coding: utf-8 -*-

from numpy import arctan2, array, cos, exp, sin, sqrt
from .constants import ANGVEL, AU_M, DAY_S, T0, pi, tau
from .earthlib import refract
from .framelib import itrs
from .functions import _T, mxm, mxv, rot_y, rot_z
from .descriptorlib import reify
from .units import Angle, Distance, _ltude
from .vectorlib import VectorFunction

class ITRSPosition(VectorFunction):
    """An x,y,z position in the Earth-centered Earth-fixed (ECEF) ITRS frame."""

    center = 399

    def __init__(self, itrs_xyz):
        self.itrs_xyz = itrs_xyz
        x, y, z = itrs_xyz.au
        self._velocity_au_per_d = ANGVEL * DAY_S * array((-y, x, 0.0))

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
        return r, v, r, None

class GeographicPosition(ITRSPosition):
    """The position of a latitude and longitude on Earth.

    Each instance of this class holds an x,y,z vector for a geographic
    position on (or above, or below) the Earth’s surface, in the ITRS
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
        self._R_latlon = mxm(
            rot_y(latitude.radians)[::-1],  # TODO: Why "::-1"?
            rot_z(-longitude.radians),
        )

    @reify  # not @property, so users have the option of overwriting it
    def target_name(self):
        m = self.elevation.m
        e = ' elevation {0:.0f} m'.format(m) if m else ''
        return '{0} latitude {1} N longitude {2} E{3}'.format(
            self.model.name, self.latitude, self.longitude, e)

    def lst_hours_at(self, t):
        """Return this position’s Local Sidereal Time in hours at time ``t``."""
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

class Geoid(object):
    """An Earth ellipsoid; maps latitudes and longitudes to x,y,z positions."""

    def __init__(self, name, radius_m, inverse_flattening):
        self.name = name
        self.radius = Distance(m=radius_m)
        self.inverse_flattening = inverse_flattening
        omf = (inverse_flattening - 1.0) / inverse_flattening
        self._one_minus_flattening_squared = omf * omf

    def latlon(self, latitude_degrees, longitude_degrees, elevation_m=0.0,
               cls=GeographicPosition):
        """Return the geographic position of a given latitude and longitude.

        Longitude is positive towards the east, so supply a negative
        number for west::

            from skyfield.api import wgs84
            observatory = wgs84.latlon(37.3414, -121.6429)

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

        ach = radius_au * c + elevation_au
        ash = radius_au * s + elevation_au

        ac = ach * cosphi
        acsst = ac * sin(lon)
        accst = ac * cos(lon)
        r = array((accst, acsst, ash * sinphi))

        return cls(self, latitude, longitude, elevation, Distance(au=r))

    def subpoint(self, position):
        """Return Earth latitude and longitude beneath a celestial ``position``.

        The input ``position`` should have a center of 399, the
        geocenter.  The return value is a `GeographicPosition` whose
        ``latitude`` and ``longitude`` are the spot on the Earth’s
        surface directly beneath the given ``position``, and whose
        ``elevation`` is the position’s distance above (or depth below)
        mean sea level.

        The underlying computation is based on Dr. T.S. Kelso's quite
        helpful article `Orbital Coordinate Systems, Part III
        <https://www.celestrak.com/columns/v02n03/>`_.

        """
        if position.center != 399:
            raise ValueError(
                'a geographic subpoint can only be calculated for positions'
                ' measured from 399, the center of the Earth, but this'
                ' position has center {0}'.format(position.center)
            )
        xyz_au = position.frame_xyz(itrs).au
        x, y, z = xyz_au

        R = sqrt(x*x + y*y)
        lon = (arctan2(y, x) - pi) % tau - pi
        lat = arctan2(z, R)

        a = self.radius.au
        f = 1.0 / self.inverse_flattening
        e2 = 2.0*f - f*f
        C = 1.0
        for iteration in 0,1,2:
            C = 1.0 / sqrt(1.0 - e2 * (sin(lat) ** 2.0))
            lat = arctan2(z + a * C * e2 * sin(lat), R)
        elevation_m = ((R / cos(lat)) - a * C) * AU_M

        return GeographicPosition(
            latitude=Angle(radians=lat),
            longitude=Angle(radians=lon),
            elevation=Distance(m=elevation_m),
            itrs_xyz=Distance(au=xyz_au),
            model=self,
        )

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

    @reify
    def R_lat(self):
        return rot_y(self.latitude.radians)[::-1]

    def itrf_xyz(self):
        return self.itrs_xyz
