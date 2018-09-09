"""Routines to solve for circumstances like sunrise, sunset, and moon phase."""

from numpy import arange, argmax, flatnonzero, linspace, multiply
from .constants import DAY_S

EPSILON = 0.001 / DAY_S

def sunrise(ephemeris, topos):
    """Return a function that, given a time, knows whether the sun is up."""
    sun = ephemeris['sun']
    topos_at = (ephemeris['earth'] + topos).at

    def is_sun_above_the_horizon(t):
        """Return `True` if the sun has risen by time `t`."""
        return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -0.8333

    return is_sun_above_the_horizon

def find(start_time, end_time, f, epsilon=EPSILON, step=None, num=12):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {} is later than your end_time {}'
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
            raise ValueError('{}{} is already true at your start time'
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
        raise ValueError('your start_time {} is later than your end_time {}'
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
