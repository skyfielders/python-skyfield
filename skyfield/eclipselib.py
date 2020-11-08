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
        # Calls to this from `find_maxima()` incur most of the expense
        # of finding eclipses, so we use raw ephemeris segments, which
        # (a) avoids computing velocities we won't use, and (b) avoids
        # computing the Earth's vector twice.

        jd, fr = t.whole, t.tdb_fraction
        b, velocity = earth_barycenter.compute_and_differentiate(jd, fr)
        e = earth.compute(jd, fr)
        m = moon.compute(jd, fr)
        s = sun.compute(jd, fr)

        earth_to_sun = s - b - e
        earth_to_moon = m - e

        # The aberration routine requires specific units.  (We can leave
        # the `earth_to_moon` vector unconverted because all we care
        # about in the end is its direction.)  We approximate the
        # Earth’s velocity as being that of the Earth-Moon barycenter.

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

    solar_radius_km = 696340
    moon_radius_km = 1.7371e3

    # Strict geometry would demand that `arcsin()` be applied to these
    # three values, but the angles are small enough that no eclipse
    # prediction seems to be affected.
    pi_m = ERAD / 1e3 / length_of(moon_to_earth)
    pi_s = ERAD / 1e3 / length_of(earth_to_sun)
    s_s = solar_radius_km / length_of(earth_to_sun)

    pi_1 = 0.998340 * pi_m

    sigma = angle_between(earth_to_sun, moon_to_earth)
    s_m = arcsin(moon_radius_km / length_of(moon_to_earth))

    penumbral = sigma < 1.02 * (pi_1 + pi_s + s_s) + s_m
    partial = sigma < 1.02 * (pi_1 + pi_s - s_s) + s_m
    total = sigma < 1.02 * (pi_1 + pi_s - s_s) - s_m

    t = t[penumbral]
    partial = partial[penumbral]
    total = total[penumbral]

    code = partial.astype(byte)
    code += total
    return t, code
