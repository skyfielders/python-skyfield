# -*- coding: utf-8 -*-
"""Raw transforms between coordinate frames, as NumPy matrices."""

from numpy import array
from .constants import ANGVEL, ASEC2RAD, DAY_S, tau
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
    """Build the matrix to rotate an ICRF vector into ecliptic coordinates."""
    _, d_eps = t._nutation_angles_radians
    true_obliquity = t._mean_obliquity_radians + d_eps
    return mxm(rot_x(- true_obliquity), t.M)

class true_equator_and_equinox_of_date(object):
    """The dynamical frame of the Earth’s true equator and equinox of date.

    This reference frame combines current theories of the Earth’s
    precession and nutation (plus a small offset between the ITRS and
    J2000 systems) to produce the orientation, on any date, of the
    Earth’s celestial poles and equator, from which a right ascension
    and declination can be produced.

    Ignoring the few tenths of an arcsecond by which the continents
    wobble with respect to the Earth’s rotational pole (“polar motion”),
    the right ascension and declination computed with this reference
    frame should reflect the real path across the sky that a body will
    describe as the Earth’s rotation carries that body along a given
    line of declination.

    See `reference_frames` for how to use frames like this one.

    """
    @staticmethod
    def rotation_at(t):
        return t.M

    @staticmethod
    def rotation_and_rate_at(t):
        # The `None` is a slight lie: t.M does have rotational velocity.
        # But it's so small that we neglect it in practice.
        return t.M, None

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

    See `reference_frames` for how to use frames like this one.

    """
    @staticmethod
    def rotation_at(t):
        R = mxm(rot_z(-t.gast * tau / 24.0), t.M)
        if t.ts.polar_motion_table is not None:
            R = mxm(t.polar_motion_matrix(), R)
        return R

    @staticmethod
    def rotation_and_rate_at(t):
        R = mxm(rot_z(-t.gast * tau / 24.0), t.M)
        if t.ts.polar_motion_table is not None:
            R = mxm(t.polar_motion_matrix(), R)
        V = _itrs_angvel_matrix
        return R, V

itrs = itrs()
