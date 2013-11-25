from numpy import sqrt


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
