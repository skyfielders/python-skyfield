# -*- coding: utf-8 -*-
"""Routines to solve for circumstances like sunrise, sunset, and moon phase."""

from __future__ import division

import numpy as np
from numpy import cos, sin, sqrt, zeros_like
from .constants import pi, tau
from .framelib import ecliptic_frame
from .searchlib import find_discrete
from .nutationlib import iau2000b_radians
from .units import Angle

_SUN = 10
_MOON = 301
_MICROSECOND = 1 / 24.0 / 3600.0 / 1e6

# Not only to support historic code but also for future convenience, let
# folks import the search routine alongside the almanac routines.
find_discrete

# Simple facts.

def phase_angle(ephemeris, body, t):
    """
    .. deprecated:: 1.42
       Use the :meth:`~skyfield.positionlib.ICRF.phase_angle()` position
       method instead.
    """
    p = ephemeris['earth'].at(t).observe(ephemeris[body])
    return p.phase_angle(ephemeris['sun'])

def fraction_illuminated(ephemeris, body, t):
    """
    .. deprecated:: 1.42
       Use the :meth:`~skyfield.positionlib.ICRF.fraction_illuminated()`
       position method instead.
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
        _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
        return (slon.radians // (tau / 4) % 4).astype(int)

    season_at.step_days = 90.0
    return season_at

MOON_PHASES = [
    'New Moon',
    'First Quarter',
    'Full Moon',
    'Last Quarter',
]

def moon_phase(ephemeris, t):
    """Return the Moon phase 0°–360° at time ``t``, where 180° is Full Moon.

    More precisely: this returns an :class:`~skyfield.units.Angle`
    giving the difference between the geocentric apparent ecliptic
    longitudes of the Moon and Sun, constrained to the interval 0°–360°
    (0–𝜏 radians) where 0° is New Moon and 180° is Full Moon.

    """
    e = ephemeris['earth'].at(t)
    moon, sun = ephemeris['moon'], ephemeris['sun']
    _, mlon, _ = e.observe(moon).apparent().frame_latlon(ecliptic_frame)
    _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
    return Angle(radians=(mlon.radians - slon.radians) % tau)

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
        _, mlon, _ = e.observe(moon).apparent().frame_latlon(ecliptic_frame)
        _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
        return ((mlon.radians - slon.radians) // (tau / 4) % 4).astype(int)

    moon_phase_at.step_days = 7.0  # one lunar phase per week
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
        lat, _, _ = e.observe(moon).apparent().frame_latlon(ecliptic_frame)
        return lat.radians > 0.0

    moon_node_at.step_days = 12.0  # 2000-2050: closest nodes 12.38 days apart
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
        _, slon, _ = e.observe(sun).apparent().frame_latlon(ecliptic_frame)
        _, tlon, _ = e.observe(target).apparent().frame_latlon(ecliptic_frame)
        return ((slon.radians - tlon.radians) / pi % 2.0).astype('int8')

    if target.target == 301:
        leading_or_trailing.step_days = 14  # Moon
    else:
        leading_or_trailing.step_days = 40  # Mercury (the fastest planet)
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

    west_of_meridian_at.step_days = 0.4  # twice a day
    return west_of_meridian_at

def sunrise_sunset(ephemeris, topos):
    """Build a function of time that returns whether the Sun is up.

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

    is_sun_up_at.step_days = 0.04  # catch days at least an hour long
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
        """Return whether the Sun is up, down, or whether there is twilight."""
        t._nutation_angles_radians = iau2000b_radians(t)
        degrees = topos_at(t).observe(sun).apparent().altaz()[0].degrees
        r = zeros_like(degrees, int)
        r[degrees >= -18.0] = 1
        r[degrees >= -12.0] = 2
        r[degrees >= -6.0] = 3
        r[degrees >= -0.8333] = 4
        return r

    is_it_dark_twilight_day_at.step_days = 0.04  # catch days at least an hour long
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

    is_body_up_at.step_days = 0.25
    return is_body_up_at

# Direct-search routines using geometry, that don't need find_discrete().

def _fastify(t):
    t._nutation_angles_radians = iau2000b_radians(t)

def _setting_hour_angle(latitude, declination, altitude_radians):
    """Return the hour angle, in radians, when a body reaches the horizon.

    Given the latitude of an observer, and the declination of a target,
    return the positive hour angle at which the body will set below the
    horizon, where the horizon is specified as `altitude_radians` above
    (positive) or below (negative) the great circle of zero altitude.

    """
    lat = latitude.radians
    dec = declination.radians
    numerator = sin(altitude_radians) - sin(lat) * sin(dec)
    denominator = cos(lat) * cos(dec)
    ha = np.arccos(np.clip(numerator / denominator, -1.0, 1.0))
    return ha

def _rising_hour_angle(latitude, declination, altitude_radians):
    return - _setting_hour_angle(latitude, declination, altitude_radians)

def _transit_ha(latitude, declination, altitude_radians):
    return 0.0

def _q(a, b, c, sign):
    discriminant = np.maximum(b*b - 4*a*c, 0.0)  # avoid tiny negative results
    return - 2*c / (b + sign * sqrt(discriminant))

def _intersection(y0, y1, v0, v1):
    # Return x at which a curve reaches y=0, given its position and
    # velocity y0,v0 at x=0 and y1,v1 at x=1.  For details, see
    # `design/intersect_function.py` in the Skyfield repository.
    sign = 1 - 2 * (y0 > y1)
    return _q(y1 - y0 - v0, v0, y0, sign)

# Per https://aa.usno.navy.mil/faq/RST_defs we estimate 34 arcminutes of
# atmospheric refraction and 16 arcminutes for the radius of the Sun.
_sun_horizon_radians = -50.0 / 21600.0 * tau
_refraction_radians = -34.0 / 21600.0 * tau
_moon_radius_m = 1.7374e6
_clip_lower = -1.0
_clip_upper = +2.0

def build_horizon_function(target):
    """Build and return a horizon function `h()` for the given `target`.

    The returned function takes a Distance argument giving the distance
    from the observer to the target, and returns a negative angle in
    radians giving the altitude which, when reached by the target's
    center, places it at the moment of rising.

    """
    target_id = getattr(target, 'target', None)
    if target_id == _SUN:
        def h(distance): return _sun_horizon_radians
    elif target_id == _MOON:
        def h(distance):
            return _refraction_radians - _moon_radius_m / distance.m
    else:
        def h(distance): return _refraction_radians
    return h

def _find(observer, target, start_time, end_time, horizon_degrees, f):
    if horizon_degrees is None:
        h = build_horizon_function(target)
    else:
        horizon_radians = horizon_degrees / 360.0 * tau
        def h(distance): return horizon_radians

    geo = observer.vector_functions[-1]  # should we check observer.center?
    latitude = geo.latitude

    # Build an array of times 0.8 days apart, in the hopes that nothing
    # ever rises (or sets or transits) twice within a 0.8-day period.
    ts = start_time.ts
    tt0 = start_time.tt
    tt1 = end_time.tt
    sample_count = int(np.ceil((tt1 - tt0) / 0.8)) + 1
    t = ts.tt_jd(np.linspace(tt0, tt1, sample_count))

    # Determine the target's hour angle and declination at those times.
    _fastify(t)
    ha, dec, distance = observer.at(t).observe(target).apparent().hadec()

    # Invoke our geometry formula: for each time `t`, predict the hour
    # angle at which the target will next reach the horizon, if its
    # declination were to remain constant.
    desired_ha_radians = f(latitude, dec, h(distance))

    # So at each time `t`, how many radians must the sky turn to bring
    # the target to the horizon?
    difference = desired_ha_radians - ha.radians
    difference %= tau

    # We want to return each rising exactly once, so where there are
    # runs of several times `t` that all precede the same rising, let's
    # throw the first few out and keep only the last one.
    i, = np.nonzero(np.diff(difference) > 0.0)

    # Trim a few arrays down to just the matching elements.
    old_ha_radians = ha.radians[i]
    old_t = t[i]

    # When might each rising have actually taken place?  Let's
    # interpolate between the two times that bracket each rising.
    a = difference[i]
    b = tau - difference[i + 1]
    tt = t.tt
    interpolated_tt = (b * tt[i] + a * tt[i+1]) / (a + b)
    t = ts.tt_jd(interpolated_tt)

    def normalize_zero_to_tau(radians):
        return radians % tau

    def normalize_plus_or_minus_pi(radians):
        return (radians + pi) % tau - pi

    normalize = normalize_zero_to_tau

    # How did we decide on this many iterations?  We played with the
    # script ./design/test_sunrise_moonrise.py in the repository.
    for i in 0, 1, 2:
        _fastify(t)

        # Expensive: generate true ha/dec at `t`.
        apparent = observer.at(t).observe(target).apparent()
        ha, dec, distance = apparent.hadec()

        # Estimate where the horizon-crossing is.
        desired_ha = f(latitude, dec, h(distance))
        ha_adjustment = desired_ha - ha.radians
        ha_adjustment = (ha_adjustment + pi) % tau - pi

        # Figure out how fast the target's HA is changing.  After two
        # iterations, just keep using the same value, in case t_diff
        # gets so small that ha_diff drops to zero.
        if i < 2:
            ha_diff = normalize(ha.radians - old_ha_radians)
            t_diff = t - old_t
            ha_per_day = ha_diff / t_diff

        # Remember this iteration's HA and `t` for the next iteration.
        old_ha_radians = ha.radians
        old_t = t

        # The big moment!  Carefully adjust `t` towards intersection.
        timebump = ha_adjustment / ha_per_day
        timebump[timebump == 0.0] = _MICROSECOND   # avoid divide-by-zero
        previous_t = t
        t = ts.tt_jd(t.whole, t.tt_fraction + timebump)

        # Different `normalize` for all but the first iteration.
        normalize = normalize_plus_or_minus_pi

    if f is _transit_ha:
        return t

    # In almost all cases, we are now happy.  But for rising and setting
    # calculations at high latitudes where the target barely scrapes the
    # horizon, we might be stuck between two solutions, and need to
    # interpolate between them.

    # Snag the observer's GeographicPosition and learn the target's
    # altitude vs the horizon and how fast it's moving vertically at
    # the second-to-last `t` we computed above.
    v = observer.vector_functions[-1]
    altitude0, _, distance0, rate0, _, _ = (
        apparent.frame_latlon_and_rates(v))

    # Even faster than _fastify(t) is to just assume that nutation
    # doesn't have much time to move over this short interval.
    t.M = previous_t.M
    t._nutation_angles_radians = previous_t._nutation_angles_radians

    # And again, this time with the very final `t` we computed.
    apparent = observer.at(t).observe(target).apparent()
    altitude1, _, distance1, rate1, _, _ = (
        apparent.frame_latlon_and_rates(v))

    # Using the target's altitude and altitude-velocity at the final
    # two times we computed, compute where it crosses the horizon.
    tdiff = t - previous_t
    t_scaled_offset = _intersection(
        altitude0.radians - h(distance0),
        altitude1.radians - h(distance1),
        rate0.radians.per_day * tdiff,
        rate1.radians.per_day * tdiff,
    )

    # In case the parabola for some reason goes crazy, don't let our
    # solution be thrown too far away from our final two times.
    t_scaled_offset = np.clip(t_scaled_offset, _clip_lower, _clip_upper)

    t = previous_t + t_scaled_offset * tdiff

    is_above_horizon =  (
        (desired_ha % pi != 0.0)
        | ((t_scaled_offset > _clip_lower) & (t_scaled_offset < _clip_upper))
    )

    return t, is_above_horizon

def find_risings(observer, target, start_time, end_time, horizon_degrees=None):
    """Return the times at which a target rises above the eastern horizon.

    Given an observer on the Earth’s surface, a target like the Sun or
    Moon or a planet, and start and stop :class:`~skyfield.timelib.Time`
    objects, this returns two arrays that have the same length.  The
    first is a :class:`~skyfield.timelib.Time` listing the moments at
    which the target rises.  The second array has ``True`` for each time
    the target really crosses the horizon, and ``False`` when the target
    merely transits without actually touching the horizon.

    See `risings-and-settings` for examples, and `horizon_degrees` for
    how to use the ``horizon_degrees`` argument.

    .. versionadded:: 1.47

    """
    return _find(observer, target, start_time, end_time, horizon_degrees,
                 _rising_hour_angle)

def find_settings(observer, target, start_time, end_time, horizon_degrees=None):
    """Return the times at which a target sets below the western horizon.

    Given an observer on the Earth’s surface, a target like the Sun or
    Moon or a planet, and start and stop :class:`~skyfield.timelib.Time`
    objects, this returns two arrays that have the same length.  The
    first is a :class:`~skyfield.timelib.Time` listing the moments at
    which the target sets.  The second array has ``True`` for each time
    the target really crosses the horizon, and ``False`` when the target
    merely transits without actually touching the horizon.

    See `risings-and-settings` for examples, and `horizon_degrees` for
    how to use the ``horizon_degrees`` argument.

    .. versionadded:: 1.47

    """
    return _find(observer, target, start_time, end_time, horizon_degrees,
                 _setting_hour_angle)

def find_transits(observer, target, start_time, end_time):
    """Return the times at which a target transits across the meridian.

    Given an observer on the Earth’s surface, a target like the Sun or
    Moon or a planet, and start and stop :class:`~skyfield.timelib.Time`
    objects, this returns a :class:`~skyfield.timelib.Time` array
    listing the moments at which the target transits across the
    meridian.

    See `transits` for example code.

    .. versionadded:: 1.47

    """
    return _find(observer, target, start_time, end_time, 0.0, _transit_ha)
