"""Generating curve parameters."""

def build_spline_given_ends(lower, upper, y0, slope0, y1, slope1):
    width = upper - lower
    slope0 = slope0 * width
    slope1 = slope1 * width
    a0 = y0
    a1 = slope0
    a2 = -2*slope0 - slope1 - 3*y0 + 3*y1
    a3 = slope0 + slope1 + 2*y0 - 2*y1
    return lower, upper, a3, a2, a1, a0
