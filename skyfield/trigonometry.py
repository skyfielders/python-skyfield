"""Routines whose currency is Angle objects."""

from numpy import arctan2, sin, cos, tan
from skyfield.constants import tau
from skyfield.units import Angle

def position_angle_of(latlon1, latlon2):
    """Return the position angle of one position with respect to another.

    Both arguments should be a tuple of :class:`~skyfield.units.Angle`
    objects whose first item is latitude-like and whose second item is
    longitude-like.  The position angle returned will be 0 degrees if
    position #2 is directly north of position #1 on the celestial
    sphere, 90 degrees if east, 180 if south, and 270 if west.

    >>> from skyfield.trigonometry import position_angle_of
    >>> from skyfield.units import Angle
    >>> a = Angle(degrees=0), Angle(degrees=0)
    >>> b = Angle(degrees=1), Angle(degrees=1)
    >>> position_angle_of(a, b)
    <Angle 315deg 00' 15.7">

    The tuples provided as arguments can have more than 2 items; any
    extra items are ignored.  This means you can pass in the output of
    routines like :meth:`~skyfield.positionlib.ICRF.ecliptic_latlon()`
    and :meth:`~skyfield.positionlib.Apparent.altaz()` as arguments.
    Their angles will be used, but the third item of each tuple, the
    distance, will be cleanly ignored.

    """
    lat1 = latlon1[0].radians
    lon1 = latlon1[1].radians
    lat2 = latlon2[0].radians
    lon2 = latlon2[1].radians
    return Angle(radians=arctan2(
        sin(lon1 - lon2),
        cos(lat1) * tan(lat2) - sin(lat1) * cos(lon2 - lon1),
    ) % tau)
