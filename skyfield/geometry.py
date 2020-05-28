"""Routines solving basic geometric problems in astronomy."""

from numpy import nan, where
from .functions import length_of

def intersect_line_and_sphere(endpoint, center, radius):
    """Compute distance to intersections of a line and a sphere.

    Given a line through the origin (0,0,0) and an (x,y,z) ``endpoint``,
    and a sphere with the (x,y,z) ``center`` and scalar ``radius``,
    return the distance from the origin to their two intersections.

    If the line is tangent to the sphere, the two intersections will be
    at the same distance.  If the line does not intersect the sphere,
    two ``nan`` values will be returned.

    """
    # See http://paulbourke.net/geometry/circlesphere/index.html#linesphere
    # Names "b" and "c" designate the familiar values from the quadratic
    # formula; happily, a = 1 because we use a unit vector for the line.

    minus_b = 2.0 * (endpoint / length_of(endpoint) * center).sum(axis=0)
    c = (center * center).sum(axis=0) - radius * radius
    discriminant = minus_b * minus_b - 4 * c
    dsqrt = discriminant ** where(discriminant < 0, nan, 0.5)  # avoid sqrt(<0)
    return (minus_b - dsqrt) / 2.0, (minus_b + dsqrt) / 2.0
