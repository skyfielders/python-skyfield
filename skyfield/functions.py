from numpy import arcsin, arctan2, array, cos, ones_like, sin, sqrt, zeros_like
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

def to_polar(xyz):
    """Convert ``[x y z]`` into spherical coordinates ``(r, theta, phi)``.

    ``r`` - vector length
    ``theta`` - angle above (+) or below (-) the xy-plane
    ``phi`` - angle around the z-axis

    The order of the three return values is intended to match ISO 31-11.

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

    The order of the three arguments is intended to match ISO 31-11.

    """
    rxy = r * cos(phi)
    return array((rxy * cos(theta), rxy * sin(theta), r * sin(phi)))

def spin_x(theta):
    z = zeros_like(theta)
    u = ones_like(theta)
    c = cos(theta)
    s = sin(theta)
    return array(((c, -s, z), (s, c, z), (z, z, u)))

def rot_x(theta):
    c = cos(theta)
    s = sin(theta)
    return array([(1.0, 0.0, 0.0), (0.0, c, s), (0.0, -s, c)])

def rot_y(theta):
    c = cos(theta)
    s = sin(theta)
    return array([(c, 0.0, -s), (0.0, 1.0, 0.0), (s, 0.0, c)])

def rot_z(theta):
    c = cos(theta)
    s = sin(theta)
    zero = theta * 0.0
    one = zero + 1.0
    return array(((c, -s, zero), (s, c, zero), (zero, zero, one)))
