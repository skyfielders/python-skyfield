# -*- coding: utf-8 -*-
"""Raw transforms between coordinate frames, as NumPy matrices."""

from numpy import array
from .constants import ANGVEL, ASEC2RAD, DAY_S, tau
from .data.spice import inertial_frames as _inertial_frames
from .functions import mxm, rot_x, rot_z

def build_matrix():
    # 'xi0', 'eta0', and 'da0' are ICRS frame biases in arcseconds taken
    # from IERS (2003) Conventions, Chapter 5.

    xi0  = -0.0166170 * ASEC2RAD
    eta0 = -0.0068192 * ASEC2RAD
    da0  = -0.01460   * ASEC2RAD

    # Compute elements of rotation matrix.

    yx = -da0
    zx =  xi0
    xy =  da0
    zy =  eta0
    xz = -xi0
    yz = -eta0

    # Include second-order corrections to diagonal elements.

    xx = 1.0 - 0.5 * (yx * yx + zx * zx)
    yy = 1.0 - 0.5 * (yx * yx + zy * zy)
    zz = 1.0 - 0.5 * (zy * zy + zx * zx)

    return array(((xx, xy, xz), (yx, yy, yz), (zx, zy, zz)))

ICRS_to_J2000 = build_matrix()
del build_matrix

def build_ecliptic_matrix(t):
    # Build the matrix to rotate an ICRF vector into ecliptic coordinates.
    _, d_eps = t._nutation_angles_radians
    true_obliquity = t._mean_obliquity_radians + d_eps
    return mxm(rot_x(- true_obliquity), t.M)

class true_equator_and_equinox_of_date(object):
    """The dynamical frame of the Earth’s true equator and equinox of date.

    This is supplied as an explicit reference frame in case you want
    x,y,z coordinates; if you want angles, it’s better to use the
    standard position method ``radec(epoch='date')`` since that will
    return the conventional units of hours-of-right-ascension instead of
    the degrees-of-longitude that ``frame_latlon()`` would return.

    This reference frame combines current theories of the Earth’s
    precession and nutation with a small offset between the ITRS and
    J2000 systems to produce right ascension and declination for a given
    date relative to the Earth’s axis and equator of rotation.

    """
    @staticmethod
    def rotation_at(t):
        return t.M

true_equator_and_equinox_of_date = true_equator_and_equinox_of_date()

_itrs_angvel_matrix = array((
    (0.0, DAY_S * ANGVEL, 0.0),
    (-DAY_S * ANGVEL, 0.0, 0.0),
    (0.0, 0.0, 0.0),
))

class itrs(object):
    """The International Terrestrial Reference System (ITRS).

    This is the IAU standard for an Earth-centered Earth-fixed (ECEF)
    coordinate system, anchored to the Earth’s crust and continents.
    This reference frame combines three other reference frames: the
    Earth’s true equator and equinox of date, the Earth’s rotation with
    respect to the stars, and (if your ``Timescale`` has polar offsets
    loaded) the polar wobble of the crust with respect to the Earth’s
    pole of rotation.

    """
    @staticmethod
    def rotation_at(t):
        R = mxm(rot_z(-t.gast * tau / 24.0), t.M)
        if t.ts.polar_motion_table is not None:
            R = mxm(t.polar_motion_matrix(), R)
        return R

    @staticmethod
    def _dRdt_times_RT_at(t):
        # TODO: taking the derivative of the instantaneous angular
        # velocity provides a more accurate transform.
        return _itrs_angvel_matrix

itrs = itrs()

class ecliptic_frame(object):
    """Reference frame of the true ecliptic and equinox of date."""
    def rotation_at(self, t):
        return build_ecliptic_matrix(t)

ecliptic_frame = ecliptic_frame()

class InertialFrame(object):
    def __init__(self, doc, matrix):
        self.__doc__ = doc
        self._matrix = matrix

    def rotation_at(self, t):
        return self._matrix

ecliptic_J2000_frame = InertialFrame(
    'Reference frame of the true ecliptic and equinox at J2000.',
    _inertial_frames['ECLIPJ2000'],
)
galactic_frame = InertialFrame(
    'Galactic System II reference frame.',
    _inertial_frames['GALACTIC'],
)
