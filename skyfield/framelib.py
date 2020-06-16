"""Raw transforms between coordinate frames, as NumPy matrices."""

from numpy import array
from .constants import ASEC2RAD
from .functions import mxm, rot_x

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
