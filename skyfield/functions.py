"""Basic operations that are needed repeatedly throughout Skyfield."""

from numpy import (
    arcsin, arctan2, array, cos, einsum, full_like,
    float_, load, rollaxis, sin, sqrt,
)
from pkgutil import get_data
from skyfield.constants import tau

def dots(v, u):
    """Given one or more vectors in `v` and `u`, return their dot products.

    This works whether `v` and `u` each have the shape ``(3,)``, or
    whether they are each whole arrays of corresponding x, y, and z
    coordinates and have shape ``(3, N)``.

    """
    return (v * u).sum(axis=0)

def _T(M):
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

_mxv, _mxm, _mxmxm = mxv, mxm, mxmxm  # In case anyone imported old name

def length_of(xyz):
    """Given a 3-element array ``[x y z]``, return its length.

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
    """Convert ``[x y z]`` to spherical coordinates ``(r, theta, phi)``.

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    Note that ``theta`` is an elevation angle measured up and down from
    the xy-plane, not a polar angle measured from the z-axis, to match
    the convention for both latitude and declination.

    """
    r = length_of(xyz)
    x, y, z = xyz
    theta = arcsin(z / r)
    phi = arctan2(y, x) % tau
    return r, theta, phi

def from_spherical(r, theta, phi):
    """Convert ``(r, theta, phi)`` to Cartesian coordinates ``[x y z]``.

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

def _to_array(value):
    """Convert plain Python sequences into NumPy arrays.

    This helps Skyfield endpoints convert caller-provided tuples and
    lists into NumPy arrays.  If the ``value`` is not a sequence, then
    it is coerced to a Numpy float object, but not an actual array.

    """
    if hasattr(value, 'shape'):
        return value
    elif hasattr(value, '__len__'):
        return array(value)
    else:
        return float_(value)

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
