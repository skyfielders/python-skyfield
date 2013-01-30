"""Raw transforms between coordinate frames, as NumPy matrices."""

from numpy import array

ASEC2RAD = 4.848136811095359935899141e-6

# 'xi0', 'eta0', and 'da0' are ICRS frame biases in arcseconds taken
# from IERS (2003) Conventions, Chapter 5.

xi0  = -0.0166170
eta0 = -0.0068192
da0  = -0.01460

# Compute elements of rotation matrix to first order the first time this
# function is called.  Elements will be saved for future use and not
# recomputed.

xx =  1.0
yx = -da0  * ASEC2RAD
zx =  xi0  * ASEC2RAD
xy =  da0  * ASEC2RAD
yy =  1.0
zy =  eta0 * ASEC2RAD
xz = -xi0  * ASEC2RAD
yz = -eta0 * ASEC2RAD
zz =  1.0

# Include second-order corrections to diagonal elements.

xx = 1.0 - 0.5 * (yx * yx + zx * zx)
yy = 1.0 - 0.5 * (yx * yx + zy * zy)
zz = 1.0 - 0.5 * (zy * zy + zx * zx)

# The actual matrices.

J2000_to_ICRS = array(((xx, xy, xz), (yx, yy, yz), (zx, zy, zz)))
ICRS_to_J2000 = array(((xx, yx, zx), (xy, yy, zy), (xz, yz, zz)))
