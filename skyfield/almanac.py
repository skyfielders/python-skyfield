"""Routines to solve for circumstances like sunrise, sunset, and moon phase."""

from numpy import cos, diff, flatnonzero, linspace, multiply, sign
from .constants import DAY_S, tau
from .nutationlib import iau2000b

EPSILON = 0.001 / DAY_S

# Simple facts.

def phase_angle(ephemeris, body, t):
    earth = ephemeris['earth']
    sun = ephemeris['sun']
    body = ephemeris[body]
    pe = earth.at(t).observe(body)
    pe.position.au *= -1     # rotate 180 degrees to point back at Earth
    t2 = t.ts.tt_jd(t.tt - pe.light_time)
    ps = body.at(t2).observe(sun)
    return pe.separation_from(ps)

def fraction_illuminated(ephemeris, body, t):
    a = phase_angle(ephemeris, body, t).radians
    return 0.5 * (1.0 + cos(a))

# Search routines.

def find_discrete(start_time, end_time, f, epsilon=EPSILON, num=12):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {0} is later than your end_time {1}'
                         .format(start_time, end_time))

    jd = linspace(jd0, jd1, (jd1 - jd0) / f.rough_period * num // 1.0)

    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        y = f(t)

        indices = flatnonzero(diff(y))
        if not len(indices):
            raise ValueError('cannot find a change in that range')

        starts = jd.take(indices)
        ends = jd.take(indices + 1)

        # Since we start with equal intervals, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(ends), y.take(indices)

def _find_maxima(start_time, end_time, f, epsilon=EPSILON, num=12):
    ts = start_time.ts
    jd0 = start_time.tt
    jd1 = end_time.tt
    if jd0 >= jd1:
        raise ValueError('your start_time {0} is later than your end_time {1}'
                         .format(start_time, end_time))

    jd = linspace(jd0, jd1, (jd1 - jd0) / f.rough_period * num // 1.0)

    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        y = f(t)

        indices = flatnonzero(diff(sign(diff(y))) == -2)
        if not len(indices):
            raise ValueError('cannot find a maximum in that range')

        starts = jd.take(indices)
        ends = jd.take(indices + 2)

        # Since we start with equal intervals, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(ends), y.take(indices)

# Discrete circumstances to search.

def sunrise_sunset(ephemeris, topos):
    """Return a function of time that returns whether the sun is up."""
    sun = ephemeris['sun']
    topos_at = (ephemeris['earth'] + topos).at

    def is_sun_up_at(t):
        """Return `True` if the sun has risen by time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -0.8333

    is_sun_up_at.rough_period = 0.5  # twice a day
    return is_sun_up_at

MOON_QUARTER_NAMES = [
    'New Moon',
    'First Quarter',
    'Full Moon',
    'Last Quarter',
]

def moon_quarter(ephemeris):
    """Return a function of time that returns the moon phase 0 through 3."""
    earth = ephemeris['earth']
    moon = ephemeris['moon']
    sun = ephemeris['sun']

    def moon_quarter_at(t):
        """Return the quarter of the moon 0 through 3 at time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        e = earth.at(t)
        _, mlon, _ = e.observe(moon).apparent().ecliptic_latlon('date')
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return ((mlon.radians - slon.radians) // (tau / 4) % 4).astype(int)

    moon_quarter_at.rough_period = 7.0  # one lunar quarter per week
    return moon_quarter_at

def _distance_to(center, target):
    def distance_at(t):
        t._nutation_angles = iau2000b(t.tt)
        distance = center.at(t).observe(target).distance().au
        return distance
    return distance_at
