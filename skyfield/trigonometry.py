"""Routines whose currency is Angle objects."""

from numpy import arctan2, sin, cos, tan
from skyfield.constants import tau
from skyfield.units import Angle

def position_angle_of(anglepair1, anglepair2):
    """Return the position angle of one position with respect to another.

    Each argument should be a tuple whose first two items are
    :class:`~skyfield.units.Angle` objects, like the tuples returned by
    :meth:`~skyfield.positionlib.ICRF.radec()`,
    :meth:`~skyfield.positionlib.ICRF.frame_latlon()`, and
    :meth:`~skyfield.positionlib.ICRF.altaz()`.

    If one of the two angle items is signed (if its ``.signed``
    attribute is true), then that angle is used as the latitude and the
    other as the longitude; otherwise the first argument is assumed to
    be the latitude.

    If the longitude has a ``.preference`` attribute of ``'hours'``, it
    is assumed to be a right ascension which is positive towards the
    east.  The position angle returned will be 0 degrees if position #2
    is directly north of position #1 on the celestial sphere, 90 degrees
    if east, 180 if south, and 270 if west.

    Otherwise, the longitude is assumed to be azimuth, which measures
    north to east around the horizon.  The position angle returned will
    be 0 degrees if position #2 is directly above position #1 in the
    sky, 90 degrees to its left, 180 if below, and 270 if to the right.

    >>> from skyfield.trigonometry import position_angle_of
    >>> from skyfield.units import Angle
    >>> a = Angle(degrees=0), Angle(degrees=0)
    >>> b = Angle(degrees=1), Angle(degrees=1)
    >>> position_angle_of(a, b)
    <Angle 315deg 00' 15.7">

    """
    lat1, lon1 = anglepair1[0], anglepair1[1]
    if lon1.signed:
        lat1, lon1 = lon1, lat1

    lat2, lon2 = anglepair2[0], anglepair2[1]
    if lon2.signed:
        lat2, lon2 = lon2, lat2

    lat1 = lat1.radians
    lon1 = lon1.radians if lon1.preference == 'hours' else -lon1.radians
    lat2 = lat2.radians
    lon2 = lon2.radians if lon2.preference == 'hours' else -lon2.radians

    return Angle(radians=arctan2(
        sin(lon2 - lon1),
        cos(lat1) * tan(lat2) - sin(lat1) * cos(lon2 - lon1),
    ) % tau)
