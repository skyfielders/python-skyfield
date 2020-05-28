"""Routines solving basic geometric problems in astronomy."""

from numpy import nan, where
from .functions import length_of

def intersect_line_and_sphere(line_endpoint, sphere_center, sphere_radius):
    """Compute distance to intersections of a line and a sphere.

    Given a line that passes through the origin (0,0,0) and an (x,y,z)
    ``line_endpoint``, and a sphere with the (x,y,z) ``sphere_center``
    and scalar ``sphere_radius``, return the distance from the origin to
    their two intersections.  If the line is tangent to the sphere, the
    two intersections will be at the same distance.  If the line does
    not intersect the sphere, two ``nan`` values will be returned.

    """
    # See http://paulbourke.net/geometry/circlesphere/index.html#linesphere
    # Names "b" and "c" designate the familiar values from the quadratic
    # formula; happily, a = 1 because we use a unit vector for the line.

    unit_vector = line_endpoint / length_of(line_endpoint)
    minus_b = 2.0 * (unit_vector * sphere_center).sum(axis=0)
    r_squared = sphere_radius * sphere_radius
    c = (sphere_center * sphere_center).sum(axis=0) - r_squared
    discriminant = minus_b * minus_b - 4 * c
    exponent = where(discriminant < 0, nan, 0.5)
    dsqrt = discriminant ** exponent
    print(discriminant, (minus_b - dsqrt) / 2.0, (minus_b + dsqrt) / 2.0)
    return (minus_b - dsqrt) / 2.0, (minus_b + dsqrt) / 2.0
