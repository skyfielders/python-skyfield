"""Routines to solve for circumstances like sunrise, sunset, and moon phase."""

from numpy import cos, diff, flatnonzero, linspace, multiply, sign
from .constants import DAY_S, tau
from .nutationlib import iau2000b

EPSILON = 0.001 / DAY_S

# Simple facts.

def phase_angle(ephemeris, body, t):
    """Compute the phase angle of a body viewed from Earth.

    The ``body`` should be an integer or string that can be looked up in
    the given ``ephemeris``, which will also be asked to provide
    positions for the Earth and Sun.  The return value will be an
    :class:`~skyfield.units.Angle` object.

    """
    earth = ephemeris['earth']
    sun = ephemeris['sun']
    body = ephemeris[body]
    pe = earth.at(t).observe(body)
    pe.position.au *= -1     # rotate 180 degrees to point back at Earth
    t2 = t.ts.tt_jd(t.tt - pe.light_time)
    ps = body.at(t2).observe(sun)
    return pe.separation_from(ps)

def fraction_illuminated(ephemeris, body, t):
    """Compute the illuminated fraction of a body viewed from Earth.

    The ``body`` should be an integer or string that can be looked up in
    the given ``ephemeris``, which will also be asked to provide
    positions for the Earth and Sun.  The return value will be a
    floating point number between zero and one.  This simple routine
    assumes that the body is a perfectly uniform sphere.

    """
    a = phase_angle(ephemeris, body, t).radians
    return 0.5 * (1.0 + cos(a))

# Search routines.

def find_discrete(start_time, end_time, f, epsilon=EPSILON, num=12):
    """Find the times when a function changes value.

    Searches between ``start_time`` and ``end_time``, which should both
    be :class:`~skyfield.timelib.Time` objects, for the occasions where
    the function ``f`` changes from one value to another.  Use this to
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

    jd = linspace(jd0, jd1, (jd1 - jd0) / f.rough_period * num // 1.0)

    end_mask = linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    o = multiply.outer

    while True:
        t = ts.tt_jd(jd)
        y = f(t)

        indices = flatnonzero(diff(y))
        if not len(indices):
            return indices, y[0:0]

        starts = jd.take(indices)
        ends = jd.take(indices + 1)

        # Since we start with equal intervals, they all should fall
        # below epsilon at around the same time; so for efficiency we
        # only test the first pair.
        if ends[0] - starts[0] <= epsilon:
            break

        jd = o(starts, start_mask).flatten() + o(ends, end_mask).flatten()

    return ts.tt_jd(ends), y.take(indices + 1)

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

SEASONS = [
    'Spring',
    'Summer',
    'Autumn',
    'Winter',
]

SEASON_EVENTS = [
    'Vernal Equinox',
    'Summer Solstice',
    'Autumnal Equinox',
    'Winter Solstice',
]

def seasons(ephemeris):
    """Build a function of time that returns the quarter of the year.

    The function that this returns will expect a single argument that is
    a :class:`~skyfield.timelib.Time` and will return 0 through 3 for
    the seasons Spring, Summer, Autumn, and Winter.

    """
    earth = ephemeris['earth']
    sun = ephemeris['sun']

    def season_at(t):
        """Return season 0 (Spring) through 3 (Winter) at time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        e = earth.at(t)
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return (slon.radians // (tau / 4) % 4).astype(int)

    season_at.rough_period = 90.0
    return season_at

def sunrise_sunset(ephemeris, topos):

    """Build a function of time that returns whether the sun is up.

    The function that this returns will expect a single argument that is
    a :class:`~skyfield.timelib.Time` and will return ``True`` if the
    sun is up, else ``False``.

    """
    sun = ephemeris['sun']
    topos_at = (ephemeris['earth'] + topos).at

    def is_sun_up_at(t):
        """Return `True` if the sun has risen by time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -0.8333

    is_sun_up_at.rough_period = 0.5  # twice a day
    return is_sun_up_at

MOON_PHASES = [
    'New Moon',
    'First Quarter',
    'Full Moon',
    'Last Quarter',
]

def moon_phases(ephemeris):
    """Build a function of time that returns the moon phase 0 through 3.

    The function that this returns will expect a single argument that is
    a :class:`~skyfield.timelib.Time` and will return the phase of the
    moon as an integer.  See the accompanying array ``MOON_PHASES`` if
    you want to give string names to each phase.

    """
    earth = ephemeris['earth']
    moon = ephemeris['moon']
    sun = ephemeris['sun']

    def moon_phase_at(t):
        """Return the phase of the moon 0 through 3 at time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        e = earth.at(t)
        _, mlon, _ = e.observe(moon).apparent().ecliptic_latlon('date')
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return ((mlon.radians - slon.radians) // (tau / 4) % 4).astype(int)

    moon_phase_at.rough_period = 7.0  # one lunar phase per week
    return moon_phase_at

def _distance_to(center, target):
    def distance_at(t):
        t._nutation_angles = iau2000b(t.tt)
        distance = center.at(t).observe(target).distance().au
        return distance
    return distance_at
