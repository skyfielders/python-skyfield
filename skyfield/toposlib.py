from numpy import einsum

from .constants import ASEC2RAD, tau
from .earthlib import terra
from .functions import rot_x, rot_y, rot_z
from .units import Distance, Angle, _interpret_ltude
from .vectorlib import VectorFunction


class Topos(VectorFunction):
    """A vector function that knows the position of a place on Earth.

    You can specify latitude and longitude by either either building
    Skyfield :class:`~skyfield.units.Angle` objects yourself and
    supplying them to the constructor as its first two arguments, or by
    providing floating point numbers to the alternate keyword arguments
    ``latitude_degrees`` and ``longitude_degrees``.

    The ``center`` of a topos object is always ``399``, the center of
    gravity of the Earth, so every call to the ``at(t)`` method of a
    topos object returns a :class:`~skyfield.positionlib.Geocentric`
    position.

    """
    center = 399
    center_name = '399 EARTH'

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
        self.elevation = Distance(m=elevation_m)
        self.x = x
        self.y = y

        self.R_lat = rot_y(latitude.radians)[::-1]

        self.target = object()  # TODO: make this more interesting
        self.target_name = '{0} N {1} E'.format(self.latitude, self.longitude)

    def __str__(self):
        return 'Topos {0}'.format(self.target_name)

    def __repr__(self):
        return '<{0}>'.format(self)

    def _snag_observer_data(self, observer_data, t):
        observer_data.altaz_rotation = self._altaz_rotation(t)
        observer_data.elevation_m = self.elevation.m

    def _altaz_rotation(self, t):
        """Compute the rotation from the ICRF into the alt-az system."""
        R_lon = rot_z(- self.longitude.radians - t.gast * tau / 24.0)
        return einsum('ij...,jk...,kl...->il...', self.R_lat, R_lon, t.M)

    def _at(self, t):
        """Compute the GCRS position and velocity of this Topos at time `t`."""
        pos, vel = terra(self.latitude.radians, self.longitude.radians,
                         self.elevation.au, t.gast)
        pos = einsum('ij...,j...->i...', t.MT, pos)
        vel = einsum('ij...,j...->i...', t.MT, vel)
        if self.x:
            R = rot_y(self.x * ASEC2RAD)
            pos = einsum('ij...,j...->i...', R, pos)
        if self.y:
            R = rot_x(self.y * ASEC2RAD)
            pos = einsum('ij...,j...->i...', R, pos)
        # TODO: also rotate velocity

        return pos, vel, pos, None
