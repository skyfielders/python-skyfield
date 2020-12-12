# -*- coding: utf-8 -*-

from numpy import array, cos, exp, sin, sqrt
from .constants import ANGVEL, DAY_S, T0
from .earthlib import refract
from .framelib import itrs
from .functions import _T, mxm, mxv, rot_y, rot_z
from .descriptorlib import reify
from .units import Angle, Distance, _ltude
from .vectorlib import VectorFunction

class GeographicPosition(VectorFunction):
    """The position of a latitude and longitude on Earth.

    This class represents a geographic position on the Earth’s surface,
    measured as an x,y,z vector in the ITRS: the international standard
    for an Earth-fixed Earth-centered (ECEF) reference frame.  Instead
    of instantiating this class directly, Skyfield users usually supply
    a longitude and latitude to an ellipsoidal model of the Earth::

        from skyfield.api import wgs84
        topos = wgs84.latlon(37.3414, -121.6429)

    Once the object has been created, here are its attributes and
    methods:

    """
    center = 399

    def __init__(self, model, latitude, longitude, elevation, itrs_position):
        self.vector_name = model.name
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
        return 'latitude {0} N longitude {1} E{2}'.format(
            self.latitude, self.longitude, e)

    def _at(self, t):
        """Compute GCRS position and velocity at time `t`."""
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

class EarthEllipsoid(object):
    """An Earth model for turning latitudes and longitudes into positions."""

    def __init__(self, name, radius_m, inverse_flattening):
        self.name = name
        self.radius = Distance(m=radius_m)
        self.inverse_flattening = inverse_flattening
        omf = (inverse_flattening - 1.0) / inverse_flattening
        self._one_minus_flattening_squared = omf * omf

    def latlon(self, latitude_degrees, longitude_degrees, elevation_m=0.0,
               cls=GeographicPosition):
        """Return a geographic position for a given latitude and longitude."""
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

grs80 = EarthEllipsoid('GRS80', 6378137.0, 298.257222101)
wgs84 = EarthEllipsoid('WGS84', 6378137.0, 298.257223563)
iers2010 = EarthEllipsoid('IERS2010', 6378136.6, 298.25642)

# Compatibility with old versions of Skyfield:

class Topos(GeographicPosition):
    """Deprecated."""
    model = EarthEllipsoid('Earth', 6378136.6, 298.25642)  # IERS2010 numbers

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
        # quite clean and clever compared to its alternatives.)
        self.model.latlon(latitude_degrees, longitude_degrees, elevation_m,
                          super(Topos, self).__init__)

    @reify
    def R_lat(self):
        return rot_y(self.latitude.radians)[::-1]

    def itrf_xyz(self):
        return self.itrs_position
