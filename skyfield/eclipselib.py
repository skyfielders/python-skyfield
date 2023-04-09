# -*- coding: utf-8 -*-
"""Search for eclipses."""

from __future__ import division

from numpy import arcsin, byte
from .constants import AU_KM, C_AUDAY, ERAD
from .functions import angle_between, length_of
from .searchlib import find_maxima
from .relativity import add_aberration

LUNAR_ECLIPSES = [
    'Penumbral',
    'Partial',
    'Total',
]

def lunar_eclipses(start_time, end_time, eph):
    """Return the lunar eclipses between ``start_time`` and ``end_time``.

    Returns a three-item tuple:

    * A :class:`~skyfield.timelib.Time` giving the dates of each eclipse.
    * An integer array of codes identifying how complete each eclipse is.
    * A dictionary of further supplementary details about each eclipse.

    This routine is adapted from the Explanatory Supplement to the
    Astronomical Almanac 11.2.3.  See `lunar-eclipses` for the details
    of how to call this function.

    """
    # Calls to the inner function `f()` from `find_maxima()` incur most
    # of the expense of this routine, so we use raw ephemeris segments.
    # This (a) avoids computing velocities we won't use (calls to `at()`
    # always compute velocity), and (b) avoids computing any segments
    # twice.
    #
    # Note that we neglect light-travel time between the Earth and Moon,
    # and also light travel time from the Sun: the Sun moves so slowly
    # that a few minutes of difference in its position does not
    # meaningfully affect our eclipse predictions.

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

        # The aberration routine requires specific units.  (We can leave
        # the `earth_to_moon` vector unconverted because we only need
        # its direction.)  We approximate the Earth’s velocity as being
        # that of the Earth-Moon barycenter.

        earth_to_sun /= AU_KM
        velocity /= AU_KM
        light_travel_time = length_of(earth_to_sun) / C_AUDAY
        add_aberration(earth_to_sun, velocity, light_travel_time)

        return angle_between(earth_to_sun, earth_to_moon)

    f.step_days = 5.0
    t, y = find_maxima(start_time, end_time, f, num=4)

    jd, fr = t.whole, t.tdb_fraction
    b = earth_barycenter.compute(jd, fr)
    e = earth.compute(jd, fr)
    m = moon.compute(jd, fr)
    s = sun.compute(jd, fr)

    earth_to_sun = s - b - e
    moon_to_earth = e - m

    solar_radius_km = 696340.0
    moon_radius_km = 1737.1

    # Strict geometry would demand that `arcsin()` be applied to these
    # three values, but the angles are small enough that no eclipse
    # prediction seems to be affected.
    pi_m = ERAD / 1e3 / length_of(moon_to_earth)
    pi_s = ERAD / 1e3 / length_of(earth_to_sun)
    s_s = solar_radius_km / length_of(earth_to_sun)

    closest_approach = angle_between(earth_to_sun, moon_to_earth)
    moon_radius = arcsin(moon_radius_km / length_of(moon_to_earth))

    # Use Danjon's method for calculating enlargement of Earth's shadow.
    # See https://eclipse.gsfc.nasa.gov/LEcat5/shadow.html
    pi_1 = 1.01 * pi_m
    penumbra_radius = pi_1 + pi_s + s_s
    umbra_radius = pi_1 + pi_s - s_s

    penumbral = closest_approach < penumbra_radius + moon_radius

    # We now know which conjunctions are eclipses!  Before proceeding,
    # narrow all of our arrays to only those incidents.
    t = t[penumbral]
    closest_approach = closest_approach[penumbral]
    moon_radius = moon_radius[penumbral]
    penumbra_radius = penumbra_radius[penumbral]
    umbra_radius = umbra_radius[penumbral]

    # Calculate fraction of the Moon's diameter covered by the Earth’s shadow.
    twice_radius = 2 * moon_radius
    umbral_magnitude = umbra_radius + moon_radius - closest_approach
    umbral_magnitude /= twice_radius
    penumbral_magnitude = penumbra_radius + moon_radius - closest_approach
    penumbral_magnitude /= twice_radius

    partial = closest_approach < umbra_radius + moon_radius
    total = closest_approach < umbra_radius - moon_radius

    code = partial.astype(byte)
    code += total

    details = {
        'closest_approach_radians': closest_approach,
        'moon_radius_radians': moon_radius,
        'penumbra_radius_radians': penumbra_radius,
        'umbra_radius_radians': umbra_radius,
        'umbral_magnitude': umbral_magnitude,
        'penumbral_magnitude': penumbral_magnitude,
    }

    return t, code, details
