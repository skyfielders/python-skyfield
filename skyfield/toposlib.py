# -*- coding: utf-8 -*-

from numpy import exp
from .constants import T0
from .earthlib import refract, terra
from .framelib import itrs
from .functions import _T, mxm, mxv, rot_y, rot_z
from .descriptorlib import reify
from .units import Angle, Distance, _interpret_ltude
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

    * The arguments ``x`` and ``y`` are ignored, and are present only
      for compatibility with earlier versions of Skyfield.

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

        p, v = terra(latitude.radians, longitude.radians, elevation.au, 0.0)
        self.itrs_position = Distance(p)
        self._velocity_au_per_d = v

        self._R_latlon = mxm(
            rot_y(self.latitude.radians)[::-1],  # TODO: Why "::-1"?
            rot_z(-self.longitude.radians),
        )

    @reify
    def R_lat(self):
        # Why did I make this a public attribute back in 2015?  Drat.
        # Now it has to stay around forever in case anyone was using it,
        # even though we now use _R_latlon instead.  At least we can
        # switch to producing it on-demand.
        return rot_y(self.latitude.radians)[::-1]

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
        R = itrs.rotation_at(t)
        R = mxm(self._R_latlon, R)
        return R

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
