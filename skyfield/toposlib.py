from numpy import einsum

from .constants import ASEC2RAD, tau
from .earthlib import terra
from .errors import raise_error_for_deprecated_time_arguments
from .functions import rot_x, rot_y, rot_z
from .positionlib import Barycentric, Geocentric
from .units import Distance, Angle, _interpret_ltude


class Topos(object):
    """A specific location on the Earth's surface.



    """
    center = 399

    # TODO(1.0): document, and add to API doc.
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

        self.target = object()  # TODO: make this more interesting

        self.latitude = latitude
        self.longitude = longitude
        self.elevation = Distance(m=elevation_m)
        self.x = x
        self.y = y

        self.R_lat = rot_y(latitude.radians)[::-1]
        self.code = self
        self.segment = _Segment(self)
        self.segments = [self.segment]
        self.ephemeris = None

    def __repr__(self):
        return '<Topos {0} N, {1} E>'.format(self.latitude, self.longitude)

    def _at(self, t):
        p, v = self.segment.icrf_vector_at(t)
        return p, v

    def _snag_observer_data(self, data, t):
        data.altaz_rotation = self._altaz_rotation(t)
        data.elevation_m = self.elevation.m

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        """Compute where this Earth location was in space on a given date."""
        tpos_au, tvel_au_per_d = self.segment.icrf_vector_at(t)
        if self.ephemeris is None:
            c = Geocentric(tpos_au, tvel_au_per_d, t)
        else:
            e = self.ephemeris['earth'].at(t)
            c = Barycentric(e.position.au + tpos_au,
                            e.velocity.au_per_d + tvel_au_per_d,
                            t)
            c._gcrs_position = tpos_au
            c._gcrs_velocity = tvel_au_per_d
        c.rGCRS = tpos_au
        c.vGCRS = tvel_au_per_d
        c.topos = self
        c.ephemeris = self.ephemeris
        c.altaz_rotation = self._altaz_rotation(t)
        return c

    def _altaz_rotation(self, t):
        """Compute the rotation from the ICRF into the alt-az system."""
        R_lon = rot_z(- self.longitude.radians - t.gast * tau / 24.0)
        return einsum('ij...,jk...,kl...->il...', self.R_lat, R_lon, t.M)


class _Segment(object):
    """Generate GCRS positions for an Earth Satellite."""

    def __init__(self, topos):
        self.center = 399
        self.target = topos

    def icrf_vector_at(self, t):
        """Return the GCRS position, velocity of this Topos at `t`."""
        topos = self.target
        pos, vel = terra(topos.latitude.radians, topos.longitude.radians,
                         topos.elevation.au, t.gast)
        pos = einsum('ij...,j...->i...', t.MT, pos)
        vel = einsum('ij...,j...->i...', t.MT, vel)
        if topos.x:
            R = rot_y(topos.x * ASEC2RAD)
            pos = einsum('ij...,j...->i...', R, pos)
        if topos.y:
            R = rot_x(topos.y * ASEC2RAD)
            pos = einsum('ij...,j...->i...', R, pos)
        # TODO: also rotate velocity
        return pos, vel
