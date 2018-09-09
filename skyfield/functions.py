"""Basic operations that are needed repeatedly throughout Skyfield."""

from numpy import arcsin, arctan2, array, cos, load, sin, sqrt
from pkgutil import get_data
from skyfield.constants import tau

def dots(v, u):
    """Given one or more vectors in `v` and `u`, return their dot products.

    This works whether `v` and `u` each have the shape ``(3,)``, or
    whether they are each whole arrays of corresponding x, y, and z
    coordinates and have shape ``(3, N)``.

    """
    return (v * u).sum(axis=0)

def length_of(xyz):
    """Given a 3-element array `[x y z]`, return its length.

    The three elements can be simple scalars, or the array can be two
    dimensions and offer three whole series of x, y, and z coordinates.

    """
    return sqrt((xyz * xyz).sum(axis=0))

def angle_between(u_vec, v_vec):
    """Given 2 vectors in `v` and `u`, return the angle separating them.

    This works whether `v` and `u` each have the shape ``(3,)``, or
    whether they are each whole arrays of corresponding x, y, and z
    coordinates and have shape ``(3, N)``. The returned angle will be
    between 0 and 180 degrees.

    This formula is from Section 12 of:
    https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf

    """
    u = length_of(u_vec)
    v = length_of(v_vec)
    num = v*u_vec - u*v_vec
    denom = v*u_vec + u*v_vec
    return 2*arctan2(length_of(num), length_of(denom))

def to_polar(xyz):
    """Convert ``[x y z]`` into spherical coordinates ``(r, theta, phi)``.

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    The meaning and order of the three return values is designed to
    match both ISO 31-11 and the traditional order used by physicists.
    Mathematicians usually define ``theta`` and ``phi`` the other way
    around, and may need to use caution when using the return values.
    See: https://en.wikipedia.org/wiki/Spherical_coordinate_system

    """
    r = length_of(xyz)
    x, y, z = xyz
    theta = arcsin(z / r)
    phi = arctan2(y, x) % tau
    return r, theta, phi

def from_polar(r, theta, phi):
    """Convert ``(r, theta, phi)`` to Cartesian coordinates ``[x y z]``.

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    The meaning and order of the three polar parameters is designed to
    match both ISO 31-11 and the traditional order used by physicists.
    Mathematicians usually define ``theta`` and ``phi`` the other way
    around, and may need to use caution when calling this function.
    See: https://en.wikipedia.org/wiki/Spherical_coordinate_system

    """
    rxy = r * cos(theta)
    return array((rxy * cos(phi), rxy * sin(phi), r * sin(theta)))

def rot_x(theta):
    c = cos(theta)
    s = sin(theta)
    zero = theta * 0.0
    one = zero + 1.0
    return array(((one, zero, zero), (zero, c, -s), (zero, s, c)))

def rot_y(theta):
    c = cos(theta)
    s = sin(theta)
    return array([(c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c)])

def rot_z(theta):
    c = cos(theta)
    s = sin(theta)
    zero = theta * 0.0
    one = zero + 1.0
    return array(((c, -s, zero), (s, c, zero), (zero, zero, one)))

try:
    from io import BytesIO
except:
    from StringIO import StringIO as BytesIO

def load_bundled_npy(filename):
    """Load a binary NumPy array file that is bundled with Skyfield."""
    data = get_data('skyfield', 'data/{0}'.format(filename))
    return load(BytesIO(data))
