"""Routines to search for maxima and zero crossings."""

from __future__ import print_function, division

from numpy import concatenate, diff, flatnonzero, linspace, multiply, sign
from .constants import DAY_S
EPSILON = 0.001 / DAY_S

def find_discrete(start_time, end_time, f, epsilon=EPSILON, num=12):
    """Find the times when a function changes value.

    Search between ``start_time`` and ``end_time``, which should both be
    :class:`~skyfield.timelib.Time` objects, for the occasions where the
    function ``f`` changes from one value to another.  Use this to
    search for events like sunrise or moon phases.

    A tuple of two arrays is returned. The first array gives the times
    at which the input function changes, and the second array specifies
    the new value of the function at each corresponding time.

    This is an expensive operation as it needs to repeatedly call the
    function to narrow down the times that it changes.  It continues
    searching until it knows each time to at least an accuracy of
    ``epsilon`` Julian days.  At each step, it creates an array of
    ``num`` new points between the lower and upper bound that it has
    established for each transition.  These two values can be changed to
    tune the behavior of the search.

    """
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {0} is later than your end_time {1}'
                         .format(start_time, end_time))

    periods = (jd1 - jd0) / f.rough_period
    if periods < 1.0:
        periods = 1.0

    jd = linspace(jd0, jd1, int(periods * num))
    return _find_discrete(ts, jd, f, epsilon, num)

# TODO: pass in `y` so it can be precomputed?

def _find_discrete(ts, jd, f, epsilon, num):
    """Algorithm core, for callers that already have a `jd` vector."""
    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        y = f(t)

        indices = flatnonzero(diff(y))
        if not len(indices):
            ends = indices  # nothing found, return empty arrays
            break

        starts = jd.take(indices)
        ends = jd.take(indices + 1)

        # Since we start with equal intervals, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            y = y.take(indices + 1)
            # Keep only the last of several zero crossings that might
            # possibly be separated by less than epsilon.
            mask = concatenate(((diff(ends) > 3.0 * epsilon), (True,)))
            ends = ends[mask]
            y = y[mask]
            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(ends), y

def _find_maxima(start_time, end_time, f, epsilon, num):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    rough_period = f.rough_period

    if jd0 >= jd1:
        raise ValueError('start_time {0} is not earlier than end_time {1}'
                         .format(start_time, end_time))

    # We find maxima by investigating every point that is higher than
    # both points next to it.  This presents a problem: if the initial
    # heights are, for example, [1.7, 1.1, 0.3, ...], there might be a
    # maximum 1.8 hidden between the first two heights, but it would not
    # meet the criteria for further investigation.  To remedy this, we
    # put an extra point out beyond each end of our range, then filter
    # our final result to remove maxima that fall outside the range.
    bump = rough_period / num
    bumps = int((jd1 - jd0) / bump) + 3
    jd = linspace(jd0 - bump, jd1 + bump, bumps)

    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        y = f(t)

        indices = flatnonzero(diff(sign(diff(y))) == -2)
        if not len(indices):
            y = y.take(indices)
            ends = indices  # nothing found, return empty arrays
            break

        starts = jd.take(indices)
        ends = jd.take(indices + 2)

        # Since we start with equal intervals, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            y = y.take(indices)

            # Filter out maxima that fell slightly outside our bounds.
            keepers = (ends >= jd0) & (ends <= jd1)
            ends = ends[keepers]
            y = y[keepers]

            # Keep only the first of several maxima that are separated
            # by less than epsilon.
            if len(ends):
                mask = concatenate(((True,), diff(ends) > epsilon))
                ends = ends[mask]
                y = y[mask]

            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(ends), y


