from numpy import sqrt


def dots(v, u):
    """Given one or more vectors in `v` and `u`, return their dot products.

    This works whether `v` and `u` each have the shape ``(3)``, or
    whether they are each whole arrays of corresponding x, y, and z
    coordinates and have shape ``(3, N)``.

    """
    return (v * u).sum(axis=0)


def length(xyz):
    """Given a vector `xyz` as a 3-value array, return its length.

    The three elements of the array can either be scalars, or can
    themselves be arrays recording whose series of x, y, and z.

    """
    return sqrt((xyz * xyz).sum(axis=0))
