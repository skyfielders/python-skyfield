# -*- coding: utf-8 -*-
"""Routines to solve for circumstances like sunrise, sunset, and moon phase."""

from __future__ import print_function, division

from numpy import arcsin, cos, zeros_like, byte
from .constants import ERAD, pi, tau
from .functions import angle_between, length_of
from .searchlib import find_discrete, find_maxima
from .nutationlib import iau2000b_radians

# Not only to support historic code but also for future convenience, let
# folks import the search routine alongside the almanac routines.
find_discrete

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

SEASON_EVENTS_NEUTRAL = [
    'March Equinox',
    'June Solstice',
    'September Equinox',
    'December Solstice',
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
        t._nutation_angles_radians = iau2000b_radians(t)
        e = earth.at(t)
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return (slon.radians // (tau / 4) % 4).astype(int)

    season_at.rough_period = 90.0
    return season_at

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
        t._nutation_angles_radians = iau2000b_radians(t)
        e = earth.at(t)
        _, mlon, _ = e.observe(moon).apparent().ecliptic_latlon('date')
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return ((mlon.radians - slon.radians) // (tau / 4) % 4).astype(int)

    moon_phase_at.rough_period = 7.0  # one lunar phase per week
    return moon_phase_at

MOON_NODES = [
    'descending',
    'ascending',
]

def moon_nodes(ephemeris):
    """Build a function of time that identifies lunar nodes.

    This returns a function taking a :class:`~skyfield.timelib.Time` and
    returning ``True`` if the Moon is above the ecliptic else ``False``.
    See :ref:`lunar-nodes` for how to use this routine.

    """
    earth = ephemeris['earth']
    moon = ephemeris['moon']

    def moon_node_at(t):
        """Return the phase of the moon 0 through 3 at time `t`."""
        e = earth.at(t)
        lat, _, _ = e.observe(moon).apparent().ecliptic_latlon('date')
        return lat.radians > 0.0

    moon_node_at.rough_period = 14.0  # one node each half lunar month
    return moon_node_at

CONJUNCTIONS = [
    'conjunction',
    'opposition',
]

def oppositions_conjunctions(ephemeris, target):
    """Build a function to find oppositions and conjunctions with the Sun.

    See :ref:`oppositions-conjunctions` for how to call this routine and
    interpret the results.

    """
    earth_at = ephemeris['earth'].at
    sun = ephemeris['sun']

    def leading_or_trailing(t):
        """Return whether the target is east or west of the Sun."""
        e = earth_at(t)
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon()
        _, tlon, _ = e.observe(target).apparent().ecliptic_latlon()
        return ((slon.radians - tlon.radians) / pi % 2.0).astype('int8')

    leading_or_trailing.rough_period = 60  # Mercury
    return leading_or_trailing

MERIDIAN_TRANSITS = ['Antimeridian transit', 'Meridian transit']

def meridian_transits(ephemeris, target, topos):
    """Build a function of time for finding when a body transits the meridian.

    The returned function accepts a :class:`~skyfield.timelib.Time`
    argument and returns ``True`` if the ``target`` body is west of the
    observer’s meridian at that time, and otherwise returns ``False.``
    See :ref:`transits` for how to use this to search for a body’s
    meridian transits and antimeridian transits.

    """
    topos_at = (ephemeris['earth'] + topos).at

    def west_of_meridian_at(t):
        """Return `True` if the target is west of the observer’s meridian."""
        t._nutation_angles_radians = iau2000b_radians(t)
        # TODO: should we work to avoid computing Topos position twice?
        # We could grab its hidden GCRS vector and do the trig ourselves.
        # Or there might be something clever we can do with the two raw
        # vectors, skipping the cost of computing spherical coordinates.
        ra1, _, _ = topos.at(t).radec(epoch='date')
        ra2, _, _ = topos_at(t).observe(target).apparent().radec(epoch='date')
        return (ra1.radians - ra2.radians) % tau < pi

    west_of_meridian_at.rough_period = 0.5  # twice a day
    return west_of_meridian_at

def sunrise_sunset(ephemeris, topos):
    """Build a function of time that returns whether the sun is up.

    The function that is returned will expect a single argument that is
    a :class:`~skyfield.timelib.Time`, and will return ``True`` if the
    sun is up, else ``False``.

    Skyfield uses the same definition as the United States Naval
    Observatory: the Sun is up when its center is 0.8333 degrees below
    the horizon, which accounts for both its apparent radius of around
    16 arcminutes and also for the 34 arcminutes by which atmospheric
    refraction on average lifts the image of the Sun.

    If you need to provide a custom value for refraction, adjust the
    estimate of the Sun’s radius, or account for a vantage point above
    the Earth’s surface, see :ref:`risings-and-settings` to learn about
    the more versatile :func:`~skyfield.almanac.risings_and_settings()`
    routine.

    """
    sun = ephemeris['sun']
    topos_at = (ephemeris['earth'] + topos).at

    def is_sun_up_at(t):
        """Return `True` if the sun has risen by time `t`.

        The Sun has risen if its altitude above the horizon is greater
        than -0.8333 degrees.

        """
        t._nutation_angles_radians = iau2000b_radians(t)
        return topos_at(t).observe(sun).apparent().altaz()[0].degrees >= -0.8333

    is_sun_up_at.rough_period = 0.5  # twice a day
    return is_sun_up_at

TWILIGHTS = {
    0: 'Night',
    1: 'Astronomical twilight',
    2: 'Nautical twilight',
    3: 'Civil twilight',
    4: 'Day',
}

def dark_twilight_day(ephemeris, topos):
    """Build a function of time returning whether it is dark, twilight, or day.

    The function that this returns will expect a single argument that is
    a :class:`~skyfield.timelib.Time` and will return:

    | 0 — Dark of night.
    | 1 — Astronomical twilight.
    | 2 — Nautical twilight.
    | 3 — Civil twilight.
    | 4 — Sun is up.

    """
    sun = ephemeris['sun']
    topos_at = (ephemeris['earth'] + topos).at

    def is_it_dark_twilight_day_at(t):
        """Return whether the sun is up, down, or whether there is twilight."""
        t._nutation_angles_radians = iau2000b_radians(t)
        degrees = topos_at(t).observe(sun).apparent().altaz()[0].degrees
        r = zeros_like(degrees, int)
        r[degrees >= -18.0] = 1
        r[degrees >= -12.0] = 2
        r[degrees >= -6.0] = 3
        r[degrees >= -0.8333] = 4
        return r

    is_it_dark_twilight_day_at.rough_period = 0.5  # twice a day
    return is_it_dark_twilight_day_at

def risings_and_settings(ephemeris, target, topos,
                         horizon_degrees=-34.0/60.0, radius_degrees=0): #?
    """Build a function of time that returns whether a body is up.

    This returns a function taking a :class:`~skyfield.timelib.Time`
    argument returning ``True`` if the body’s altazimuth altitude angle
    plus ``radius_degrees`` is greater than ``horizon_degrees``, else
    ``False``.  See :ref:`risings-and-settings` to learn about how to
    search for risings and settings, and to see more about using the
    parameters ``horizon_degrees`` and ``radius_degrees``.

    """
    topos_at = (ephemeris['earth'] + topos).at
    h = horizon_degrees - radius_degrees

    def is_body_up_at(t):
        """Return `True` if the target has risen by time `t`."""
        t._nutation_angles_radians = iau2000b_radians(t)
        return topos_at(t).observe(target).apparent().altaz()[0].degrees > h

    is_body_up_at.rough_period = 0.5  # twice a day
    return is_body_up_at

LUNAR_ECLIPSES = [
    'Penumbral',
    'Partial',
    'Total',
]

from .relativity import add_aberration, C_AUDAY
from .constants import AU_KM, C

def lunar_eclipses(eph, start_time, end_time):
    """Return the lunar eclipses between `start_time` and `end_time`.

    Two arrays are returned: a :class:`~skyfield.timelib.Time` giving
    the dates of each eclipse, and an integer array of codes (whose
    meanings are listed in the `LUNAR_ECLIPSES` list).  Adapted from the
    Explanatory Supplement to the Astronomical Almanac 11.2.3.

    """
    sdict = dict(((s.center, s.target), s.spk_segment) for s in eph.segments)
    sun = sdict[0,10]
    earth_barycenter = sdict[0,3]
    earth = sdict[3,399]
    moon = sdict[3,301]

    def f(t):
        jd, fr = t.whole, t.tdb_fraction

        b, velocity = earth_barycenter.compute_and_differentiate(jd, fr)
        e = earth.compute(jd, fr)
        m = moon.compute(jd, fr)
        s = sun.compute(jd, fr)

        earth_to_sun = s - b - e
        earth_to_moon = m - e

        earth_to_sun /= AU_KM
        velocity /= AU_KM
        light_travel_time = length_of(earth_to_sun) / C_AUDAY

        add_aberration(earth_to_sun, velocity, light_travel_time)
        return angle_between(earth_to_sun, earth_to_moon)

    def f2(t):
        e = earth.at(t).position.au
        m = moon.at(t).position.au
        #s = sun.at(t).position.au
        #v1 = s - e
        #v1 = moon.at(t).observe(sun).position.au
        #v1 = moon.at(t).observe(sun).apparent().position.au
        #v1 = earth.at(t).observe(sun).position.au
        v1 = earth.at(t).observe(sun).apparent().position.au
        v2 = m - e
        #v2 = earth.at(t).observe(moon).position.au
        #v2 = -moon.at(t).observe(earth).position.au
        return angle_between(v1, v2)

    f.step_days = 5.0
    t, y = find_maxima(start_time, end_time, f)

    earth = eph['earth']
    moon = eph['moon']
    sun = eph['sun']

    e = earth.at(t).position.m
    s = sun.at(t).position.m
    m = moon.at(t).position.m

    #me = moon.at(t).observe(earth).apparent().position.m
    #me = moon.at(t).observe(earth).apparent().position.m
    me = e - m
    es = s - e
    #es = s - m
    #es = moon.at(t).observe(sun).apparent().position.m
    #es = earth.at(t).observe(sun).apparent().position.m

    solar_radius_m = 696340e3
    moon_radius_m = 1.7371e6

    # pi_m = arcsin(ERAD / length_of(me))
    # X =ERAD / length_of(me)
    # Y = arcsin(X)
    # print(X, X/Y)
    pi_m = ERAD / length_of(me)
    #pi_s = arcsin(ERAD / length_of(es))
    pi_s = ERAD / length_of(es)
    #s_s = arcsin(solar_radius_m / length_of(es))
    s_s = solar_radius_m / length_of(es)

    pi_1 = 0.998340 * pi_m

    sigma = angle_between(es, me)
    s_m = arcsin(moon_radius_m / length_of(me))

    penumbral = sigma < 1.02 * (pi_1 + pi_s + s_s) + s_m
    partial = sigma < 1.02 * (pi_1 + pi_s - s_s) + s_m
    total = sigma < 1.02 * (pi_1 + pi_s - s_s) - s_m

    t = t[penumbral]
    partial = partial[penumbral]
    total = total[penumbral]

    code = partial.astype(byte)
    code += total
    return t, code
