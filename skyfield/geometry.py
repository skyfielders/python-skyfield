"""Routines solving basic geometric problems in astronomy."""

from .functions import dots, length_of, nan, sqrt_nan, where

def intersect_line_and_sphere(endpoint, center, radius):
    """Compute distance to intersections of a line and a sphere.

    Given a line through the origin (0,0,0) and an |xyz| ``endpoint``,
    and a sphere with the |xyz| ``center`` and scalar ``radius``,
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
    dsqrt = sqrt_nan(discriminant)
    return (minus_b - dsqrt) / 2.0, (minus_b + dsqrt) / 2.0

def line_and_ellipsoid_intersection(line_start, line_direction, radii):
    """Return the |xyz| position where a line intersects an ellipsoid.

    All three arguments are |xyz| arrays.  The line is specified by a
    ``line_start`` endpoint and a ``line_direction`` vector.  The
    ellipsoid is centered at the origin and is specified by its three
    ``radii`` that point along the three coordinate axes.

    Returns the |xyz| point of intersection, or ``[nan nan nan]`` if the
    line does not intersect the sphere.

    """
    # Based on `surfpt.f` from the SPICE Toolkit.
    if len(getattr(line_start, 'shape', ())) > 1:
        radii = radii.reshape((3, 1))

    # Scale coordinates so the ellipsoid becomes the unit sphere.
    start = line_start / radii
    direction = line_direction / radii

    # Where does the line come closest to the sphere's center?
    closest_point = start - _vector_projection(start, direction)

    startmag = length_of(start)
    pmag = length_of(closest_point)

    is_inside_sphere = startmag < 1.0
    is_behind_us = dots(closest_point - start, direction) < 0.0

    sign = where(is_inside_sphere, +1.0, where(is_behind_us, nan, -1.0))
    half_chord_length = sqrt_nan(1.0 - pmag*pmag)
    unit_direction = direction / length_of(direction)
    intersection = closest_point + sign * half_chord_length * unit_direction
    return where(startmag == 1.0, line_start, intersection * radii)

def _vector_projection(a, b):
    return dots(a,b) / dots(b,b) * b
