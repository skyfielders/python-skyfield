# -*- coding: utf-8 -*-

from numpy import array, cos, exp, sin, sqrt
from .constants import ANGVEL, DAY_S, T0
from .earthlib import refract, terra
from .framelib import itrs
from .functions import _T, mxm, mxv, rot_y, rot_z
from .descriptorlib import reify
from .units import Angle, Distance, _interpret_ltude
from .vectorlib import VectorFunction

class EarthEllipsoid(object):
    def __init__(self, name, radius_m, inverse_flattening):
        self.name = name
        self.radius = Distance(m=radius_m)
        self.inverse_flattening = inverse_flattening
        omf = (inverse_flattening - 1.0) / inverse_flattening
        self._one_minus_flattening_squared = omf * omf

    def latlon(self, latitude_degrees, longitude_degrees, elevation_m=0.0):
        latitude = Angle(degrees=latitude_degrees)
        longitude = Angle(degrees=longitude_degrees)
        elevation = Distance(m=elevation_m)

        lat = latitude.radians
        lon = longitude.radians
        radius_au = self.radius.au
        elevation_au = elevation.au

        sinphi = sin(lat)
        cosphi = cos(lat)
        c = 1.0 / sqrt(cosphi * cosphi +
                       sinphi * sinphi * self._one_minus_flattening_squared)
        s = self._one_minus_flattening_squared * c

        ach = radius_au * c + elevation_au
        ash = radius_au * s + elevation_au

        sinst = sin(lon)
        cosst = cos(lon)

        ac = ach * cosphi
        acsst = ac * sinst
        accst = ac * cosst
        r = array((accst, acsst, ash * sinphi))

        return Topocentric(self, latitude, longitude, elevation, Distance(au=r))

grs80 = EarthEllipsoid('GRS80', 6378137.0, 298.257222101)
wgs84 = EarthEllipsoid('WGS84', 6378137.0, 298.257223563)
iers2010 = EarthEllipsoid('IERS2010', 6378136.6, 298.25642)

class Topocentric(VectorFunction):
    """A vector function that knows the position of a place on Earth.

    This class represents a geographic position on the Earth’s surface.
    Instead of instantiating this class directly, Skyfield users usually
    build an instance using an ellipsoidal model of the Earth::

        from skyfield.api import wgs84
        topos = wgs84.latlon(37.3414, -121.6429)

    Once the object has been created, here are its attributes and
    methods:

    """
    center = 399

    def __init__(self, model, latitude, longitude, elevation, itrs_position):
        self.model = model
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.itrs_position = itrs_position

        r = itrs_position.au
        x, y = r[0], r[1]
        self._velocity_au_per_d = ANGVEL * DAY_S * array((-y, x, 0.0))

        self._R_latlon = mxm(
            rot_y(latitude.radians)[::-1],  # TODO: Why "::-1"?
            rot_z(-longitude.radians),
        )

    @property
    def target(self):
        # When used as a vector function, this Earth geographic location
        # computes positions from the Earth's center to itself.  (This
        # is a property, rather than an attribute, to avoid a circular
        # reference that delays garbage collection.)
        return self

    @reify  # not @property, so users have the option of overwriting it
    def target_name(self):
        m = self.elevation.m
        e = ' elevation {0:.0f} m'.format(m) if m else ''
        return '{0} latitude {1} N longitude {2} E{3}'.format(
            self.model.name, self.latitude, self.longitude, e)

    def _at(self, t):
        """Compute the GCRS position and velocity of this Topos at time `t`."""
        r = self.itrs_position.au
        v = self._velocity_au_per_d

        RT = _T(itrs.rotation_at(t))
        r = mxv(RT, r)
        v = mxv(RT, v)
        return r, v, r, None

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
        """Compute rotation from ICRF to this location’s altazimuth system."""
        return mxm(self._R_latlon, itrs.rotation_at(t))

class Topos(Topocentric):
    model = EarthEllipsoid('Earth', 6378136.6, 298.25642)  # IERS2010 numbers

    def __init__(self, latitude=None, longitude=None, latitude_degrees=None,
                 longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):

        if latitude_degrees is not None:
            latitude = Angle(degrees=latitude_degrees)
        elif isinstance(latitude, (str, float, tuple)):
            latitude = _interpret_ltude(latitude, 'latitude', 'N', 'S')
        elif not isinstance(latitude, Angle):
            raise TypeError('please provide either latitude_degrees=<float>'
                            ' or latitude=<skyfield.units.Angle object>'
                            ' with north being positive')

        if longitude_degrees is not None:
            longitude = Angle(degrees=longitude_degrees)
        elif isinstance(longitude, (str, float, tuple)):
            longitude = _interpret_ltude(longitude, 'longitude', 'E', 'W')
        elif not isinstance(longitude, Angle):
            raise TypeError('please provide either longitude_degrees=<float>'
                            ' or longitude=<skyfield.units.Angle object>'
                            ' with east being positive')

        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation = Distance(m=elevation_m)

        p, v = terra(latitude.radians, longitude.radians, elevation.au, 0.0)
        self.itrs_position = Distance(p)
        self._velocity_au_per_d = v

        self._R_latlon = mxm(
            rot_y(self.latitude.radians)[::-1],  # TODO: Why "::-1"?
            rot_z(-self.longitude.radians),
        )

    @reify
    def R_lat(self):
        return rot_y(self.latitude.radians)[::-1]

    def itrf_xyz(self):
        return self.itrs_position
