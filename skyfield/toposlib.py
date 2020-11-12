# -*- coding: utf-8 -*-

from numpy import exp
from .constants import ASEC2RAD, T0, tau
from .earthlib import refract, terra
from .functions import mxmxm, mxv, rot_x, rot_y, rot_z
from .descriptorlib import reify
from .units import Angle, Distance, Velocity, _interpret_ltude
from .vectorlib import VectorFunction

class Topos(VectorFunction):
    """A vector function that knows the position of a place on Earth.

    The constructor needs:

    * Either an :class:`~skyfield.units.Angle` for the ``latitude`` or
      else a plain float ``latitude_degrees`` providing the angle in
      degrees.

    * Either an :class:`~skyfield.units.Angle` for the ``longitude`` or
      else a plain float ``longitude_degrees`` providing the angle in
      degrees.

    * Optionally, the ``elevation_m`` of the location, in meters above
      mean sea level on a WGS-84 globe.  If not specified, the location
      will be assumed to sit at exactly sea level.

    The ``center`` of a topos object is always ``399``, the center of
    gravity of the Earth, so every call to the ``at(t)`` method of a
    topos object returns a :class:`~skyfield.positionlib.Geocentric`
    position.

    Once the object has been created, here are its attributes and
    methods:

    """
    center = 399

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
        self.x = x
        self.y = y

        self.R_lat = rot_y(latitude.radians)[::-1]

        p, v = terra(latitude.radians, longitude.radians, elevation.au, 0.0)
        self.itrs_position = Distance(p)
        self.itrs_velocity = Velocity(v)

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
        return 'Earth latitude {0} N longitude {1} E{2}'.format(
            self.latitude, self.longitude, e)

    def _altaz_rotation(self, t):
        """Compute the rotation from the ICRF into the alt-az system."""
        R_lon = rot_z(- self.longitude.radians - t.gast * tau / 24.0)
        return mxmxm(self.R_lat, R_lon, t.M)

    def _at(self, t):
        """Compute the GCRS position and velocity of this Topos at time `t`."""
        pos = self.itrs_position.au
        vel = self.itrs_velocity.au_per_d

        if self.x or self.y:
            R = rot_x(-self.y * ASEC2RAD)
            pos = mxv(R, pos)

            R = rot_y(-self.x * ASEC2RAD)
            pos = mxv(R, pos)

            sprime = -47.0e-6 * (t.tdb - T0) / 36525.0
            tiolon = -sprime * ASEC2RAD
            pos = mxv(rot_z(-tiolon), pos)

        theta = t.gast / 24.0 * tau
        R = rot_z(theta)
        pos = mxv(R, pos)
        vel = mxv(R, vel)
        pos = mxv(t.MT, pos)
        vel = mxv(t.MT, vel)

        return pos, vel, pos, None

    def itrf_xyz(self):
        """DEPRECATED: access the ``itrs_position`` attribute instead."""
        return self.itrs_position

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
        """Compute the altazimuth rotation matrix for this location’s sky."""
        return self._altaz_rotation(t)
