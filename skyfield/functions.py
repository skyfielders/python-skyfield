"""Basic operations that are needed repeatedly throughout Skyfield."""

from numpy import (
    arctan2, array, cos, einsum, finfo, float64,
    full_like, hypot, load, nan, rollaxis, sin, sqrt, where,
)
from pkgutil import get_data
from skyfield.constants import tau

_AVOID_DIVIDE_BY_ZERO = finfo(float64).tiny

class A(object):
    """Allow literal NumPy arrays to be spelled ``A[1, 2, 3]``."""
    __getitem__ = array
A = A()

def sqrt_nan(n):
    """Return the square root of ``n``, or ``nan`` if ``n < 0``."""
    # (See design/sqrt_nan.py for a speed comparison of approaches.)
    return where(n < 0.0, nan, sqrt(abs(n)))

def dots(v, u):
    """Given one or more vectors in `v` and `u`, return their dot products.

    This works whether `v` and `u` each have the shape ``(3,)``, or
    whether they are each whole arrays of corresponding x, y, and z
    coordinates and have shape ``(3, N)``.

    """
    return (v * u).sum(axis=0)

def T(M):
    """Swap the first two dimensions of an array."""
    return rollaxis(M, 1)

def mxv(M, v):
    """Matrix times vector: multiply an NxN matrix by a vector."""
    return einsum('ij...,j...->i...', M, v)

def mxm(M1, M2):
    """Matrix times matrix: multiply two NxN matrices."""
    return einsum('ij...,jk...->ik...', M1, M2)

def mxmxm(M1, M2, M3):
    """Matrix times matrix times matrix: multiply 3 NxN matrices together."""
    return einsum('ij...,jk...,kl...->il...', M1, M2, M3)

_T, _mxv, _mxm, _mxmxm = T, mxv, mxm, mxmxm  # In case anyone imported old name

def length_of(xyz):
    """Given a 3-element array |xyz|, return its length.

    The three elements can be simple scalars, or the array can be two
    dimensions and offer three whole series of x, y, and z coordinates.

    """
    return sqrt((xyz * xyz).sum(axis=0))

def angle_between(u, v):
    """Given two vectors `v` and `u`, return the radian angle separating them.

    This works whether `v` and `u` each have the shape ``(3,)``, or
    whether they are each whole arrays of corresponding x, y, and z
    coordinates with shape ``(3, N)``. The returned angle will be
    between 0 and tau/2.

    This formula is from Section 12 of:
    https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf

    """
    a = u * length_of(v)
    b = v * length_of(u)
    return 2.0 * arctan2(length_of(a - b), length_of(a + b))

def to_spherical(xyz):
    """Convert |xyz| to spherical coordinates (r,theta,phi).

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    Note that ``theta`` is an elevation angle measured up and down from
    the xy-plane, not a polar angle measured from the z-axis, to match
    the convention for both latitude and declination.

    """
    r = length_of(xyz)
    x, y, z = xyz
    theta = arctan2(z, hypot(x, y))
    phi = arctan2(y, x) % tau
    return r, theta, phi

def _to_spherical_and_rates(r, v):
    # Convert Cartesian rate and velocity vectors to angles and rates.
    x, y, z = r
    xdot, ydot, zdot = v

    length = length_of(r)
    lat = arctan2(z, hypot(x, y));
    lon = arctan2(y, x) % tau
    range_rate = dots(r, v) / length

    x2 = x * x
    y2 = y * y
    x2_plus_y2 = x2 + y2 + _AVOID_DIVIDE_BY_ZERO
    lat_rate = (x2_plus_y2 * zdot - z * (x * xdot + y * ydot)) / (
        (x2_plus_y2 + z*z) * sqrt(x2_plus_y2))
    lon_rate = (x * ydot - xdot * y) / x2_plus_y2

    return length, lat, lon, range_rate, lat_rate, lon_rate

def from_spherical(r, theta, phi):
    """Convert (r,theta,phi) to Cartesian coordinates |xyz|.

    ``r`` - vector length
    ``theta`` - angle in radians above (+) or below (-) the xy-plane
    ``phi`` - angle in radians around the z-axis

    Note that ``theta`` is an elevation angle measured up and down from
    the xy-plane, not a polar angle measured from the z-axis, to match
    the convention for both latitude and declination.

    """
    rxy = r * cos(theta)
    return array((rxy * cos(phi), rxy * sin(phi), r * sin(theta)))

# Support users who might have imported these under their old names.
# I'm not sure why I called what are clearly spherical coordinates "polar".
to_polar = to_spherical
from_polar = from_spherical

def rot_x(theta):
    c = cos(theta)
    s = sin(theta)
    zero = theta * 0.0
    one = zero + 1.0
    return array(((one, zero, zero), (zero, c, -s), (zero, s, c)))

def rot_y(theta):
    c = cos(theta)
    s = sin(theta)
    zero = theta * 0.0
    one = zero + 1.0
    return array(((c, zero, s), (zero, one, zero), (-s, zero, c)))

def rot_z(theta):
    c = cos(theta)
    s = sin(theta)
    zero = theta * 0.0
    one = zero + 1.0
    return array(((c, -s, zero), (s, c, zero), (zero, zero, one)))

# The rotation matrices R1, R2, and R3 in _The Explanatory Supplement to
# the Astronomical Almanac_ use a left-handed rotation around the x and
# z axes.  In case anyone needs them:

def R1(theta): return rot_x(-theta)
def R2(theta): return rot_y(theta)
def R3(theta): return rot_z(-theta)

def angular_velocity_matrix(angular_velocity_vector):
    x, y, z = angular_velocity_vector
    zero = x * 0.0
    return array(((zero, -z, y), (z, zero, -x), (-y, x, zero)))

def _to_array(value):
    """Convert plain Python sequences into NumPy arrays.

    This lets users pass plain old Python lists and tuples to Skyfield,
    instead of always having to remember to build NumPy arrays.  We pass
    any kind of generic sequence to the NumPy ``array()`` constructor
    and wrap any other kind of value in a NumPy ``float64`` object.

    """
    if hasattr(value, 'shape'):
        return value
    elif hasattr(value, '__len__'):
        return array(value)
    else:
        return float64(value)

def _reconcile(a, b):
    """Coerce two NumPy generics-or-arrays to the same number of dimensions."""
    an = getattr(a, 'ndim', 0)
    bn = getattr(b, 'ndim', 0)
    difference = bn - an
    if difference > 0:
        if an:
            a.shape += (1,) * difference
        else:
            a = full_like(b, a)
    elif difference < 0:
        if bn:
            b.shape += (1,) * -difference
        else:
            b = full_like(a, b)
    return a, b

try:
    from io import BytesIO
except:
    from StringIO import StringIO as BytesIO

def load_bundled_npy(filename):
    """Load a binary NumPy array file that is bundled with Skyfield."""
    data = get_data('skyfield', 'data/{0}'.format(filename))
    return load(BytesIO(data))
