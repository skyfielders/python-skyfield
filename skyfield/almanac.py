"""Routines to solve for circumstances like sunrise, sunset, and moon phase."""

from numpy import arange, argmax, flatnonzero, linspace, multiply
from .constants import DAY_S, tau

EPSILON = 0.001 / DAY_S

def find(start_time, end_time, f, epsilon=EPSILON, step=None, num=12):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {0} is later than your end_time {1}'
                         .format(start_time, end_time))
    step = 0.1  # TODO

    jd = arange(jd0, jd1, step)
    while True:
        t = ts.tt_jd(jd)
        true_false = f(t)
        i = argmax(true_false)
        if i == 0:
            name = getattr(f, '__name__', None)
            parens = '' if name is None else '()'
            raise ValueError('{0}{1} is already true at your start time'
                             .format(name or 'your criterion', parens))
        jd0, jd1 = jd[i-1], jd[i]
        if jd1 - jd0 <= epsilon:
            break
        jd = linspace(jd0, jd1, num)

    return ts.tt_jd(jd0)

def find_all(start_time, end_time, f, epsilon=EPSILON, step=None, num=12):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {0} is later than your end_time {1}'
                         .format(start_time, end_time))
    step = 0.1  # TODO
    jd = arange(jd0, jd1, step)

    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        true_false = f(t)
        rising_edges = flatnonzero(~true_false[:-1] & true_false[1:])
        starts = jd.take(rising_edges)
        ends = jd.take(rising_edges + 1)

        # Since we create the intervals equal, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(starts)

def find_all2(start_time, end_time, f, epsilon=EPSILON, step=None, num=12):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {0} is later than your end_time {1}'
                         .format(start_time, end_time))
    step = 0.1  # TODO
    step = 7.0
    jd = arange(jd0, jd1, step)

    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        from numpy import diff
        y = f(t)
        indices = flatnonzero(diff(y))
        starts = jd.take(indices)
        ends = jd.take(indices + 1)

        # Since we create the intervals equal, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(ends)

def sunrise_sunset(ephemeris, topos):
    """Return a function of time that returns whether the sun is up."""
    sun = ephemeris['sun']
    topos_at = (ephemeris['earth'] + topos).at

    def is_sun_above_the_horizon_at(t):
        """Return `True` if the sun has risen by time `t`."""
        return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -0.8333

    return is_sun_above_the_horizon_at

# MOON_QUARTER_NAMES = {
#     'First Quarter',
# }

def moon_quarter(ephemeris):
    """Return a function of time that returns the moon phase 0 through 3."""
    earth = ephemeris['earth']
    moon = ephemeris['moon']
    sun = ephemeris['sun']

    def moon_quarter_at(t):
        """Return `True` if the sun has risen by time `t`."""
        e = earth.at(t)
        _, mlon, _ = e.observe(moon).apparent().ecliptic_latlon('date')
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return (mlon.radians - slon.radians) // (tau / 4) % 4

    return moon_quarter_at
