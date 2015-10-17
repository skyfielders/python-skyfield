from numpy import einsum

from .constants import ASEC2RAD, tau
from .earthlib import terra
from .ephemerislib import Body, Segment
from .functions import rot_x, rot_y, rot_z
from .positionlib import Barycentric, Geocentric
from .timelib import takes_julian_date
from .units import Distance, Angle, _interpret_ltude


class Topos(Body):
    """An object representing a specific location on the Earth's surface."""

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
        self.code = self
        self.segments = [Segment(399, self, self.compute)]
        self.ephemeris = None

    def __repr__(self):
        return '<Topos {0} N, {1} E>'.format(self.latitude, self.longitude)

    def compute(self, jd):
        position, velocity = self._position_and_velocity(jd)
        return position, velocity

    @takes_julian_date
    def at(self, jd):
        """Compute where this Earth location was in space on a given date."""
        tpos_au, tvel_au_per_d = self._position_and_velocity(jd)
        if self.ephemeris is None:
            c = Geocentric(tpos_au, tvel_au_per_d, jd)
        else:
            e = self.ephemeris['earth'].at(jd)
            c = Barycentric(e.position.au + tpos_au,
                            e.velocity.au_per_d + tvel_au_per_d,
                            jd)
            c.geocentric = False  # test, then get rid of this attribute
        c.rGCRS = tpos_au
        c.vGCRS = tvel_au_per_d
        c.topos = self
        c.ephemeris = self.ephemeris
        c.altaz_rotation = self._altaz_rotation(jd)
        return c

    def _position_and_velocity(self, jd):
        """Return the GCRS position, velocity of this Topos at `jd`."""
        pos, vel = terra(self.latitude.radians, self.longitude.radians,
                         self.elevation.au, jd.gast)
        pos = einsum('ij...,j...->i...', jd.MT, pos)
        vel = einsum('ij...,j...->i...', jd.MT, vel)
        if self.x:
            R = rot_y(self.x * ASEC2RAD)
            pos = einsum('ij...,j...->i...', R, pos)
        if self.y:
            R = rot_x(self.y * ASEC2RAD)
            pos = einsum('ij...,j...->i...', R, pos)
        # TODO: also rotate velocity
        return pos, vel

    def _altaz_rotation(self, jd):
        """Compute the rotation from the ICRS into the alt-az system."""
        R_lon = rot_z(- self.longitude.radians - jd.gast * tau / 24.0)
        return einsum('ij...,jk...,kl...->il...', self.R_lat, R_lon, jd.M)
