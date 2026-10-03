# -*- coding: utf-8 -*-
"""Search for eclipses."""

from __future__ import division

from numpy import arcsin, byte, arctan2, sin, pi, arccos, sqrt, dot, array, empty
from .constants import AU_KM, C_AUDAY, ERAD
from .functions import angle_between, length_of
from .searchlib import find_maxima, find_minima
from .relativity import add_aberration

LUNAR_ECLIPSES = [
    'Penumbral',
    'Partial',
    'Total',
]

SOLAR_ECLIPSES = [
    'Partial',
    'Total/Hybrid',
    'Annular'
]


def _compute_angle(t, earth_barycenter, earth, moon, sun):
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

    f = lambda x: _compute_angle(x, earth_barycenter, earth, moon, sun)
    f.step_days = 5.0

    t, y = find_maxima(start_time, end_time, f , num=4)

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

def _cross_area(d, r, R):
    """
    Computes common area of the circles with radius ``r`` and ``R``
    that are separated by distance ``d``.
    """

    return (
        r * r * arccos((d * d + r * r - R * R) / (2 * d * r))
        + R * R * arccos((d * d + R * R - r * r) / (2 * d * R))
        - 1 / 2 * sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R))
    )


def _line_sphere_intersection(unit_vec, sphere_center, sphere_radius):
    """
    Computes closest point of line intersecting the sphere.
    The line starts at origin (0, 0, 0) and has unit vector ``unit_vec``.
    If line does intersect the sphere, returns closest point of sphere
    surface to the line.
    """

    center_sq = dot(sphere_center, sphere_center)
    prod1 = dot(unit_vec, sphere_center)
    delta = prod1 * prod1 - center_sq + sphere_radius * sphere_radius
    if delta < 0:
        s_to_point = unit_vec * prod1
        unit_e_to_point = (s_to_point - sphere_center) / length_of(
            s_to_point - sphere_center
        )
        return sphere_center + unit_e_to_point * sphere_radius
    else:
        d = prod1 - sqrt(delta)
        return unit_vec * d


def _magnitude_obscurity(e_to_s, e_to_m, s_radius, m_radius, e_radius):
    """
    Calculates magnitude and obscurity of the eclipse.
    """

    s_to_m = e_to_m - e_to_s
    s_to_ecl = _line_sphere_intersection(s_to_m / length_of(s_to_m),
                                         -e_to_s, e_radius)
    m_to_ecl = s_to_ecl - s_to_m
    scale = length_of(s_to_ecl) / length_of(m_to_ecl)
    s_radius_res = s_radius / scale
    s_to_ecl_res = s_to_ecl / scale
    dist = length_of(m_to_ecl - s_to_ecl_res)
    if dist >= s_radius_res + m_radius:
        return 0.0, 0.0
    mag = m_radius / s_radius_res
    if m_radius + dist <= s_radius_res:
        obs = m_radius * m_radius / (s_radius_res * s_radius_res)
    elif s_radius_res + dist <= m_radius:
        obs = 1.0
    else:
        mag = (s_radius_res + m_radius - dist) / (2 * s_radius_res)
        cross_a = _cross_area(dist, s_radius_res, m_radius)
        obs = cross_a / (pi * s_radius_res * s_radius_res)
    return mag, obs


def solar_eclipses(start_time, end_time, eph):
    """Return the solar eclipses between ``start_time`` and ``end_time``.

    Returns a 4-item tuple:

    * A :class:`~skyfield.timelib.Time` giving the dates of each eclipse.
    * An integer array of codes identifying what type eclipse is.
    * A float indicating the magnitude of the eclipse.
    * A float indicating the obscurity of the eclipse.

    Comments regarding approximations for the method lunar_eclipses()
    are applicable for this routine as well.
    """

    sdict = dict(((s.center, s.target), s.spk_segment) for s in eph.segments)
    sun = sdict[0,10]
    earth_barycenter = sdict[0,3]
    earth = sdict[3,399]
    moon = sdict[3,301]

    f = lambda x: _compute_angle(x, earth_barycenter, earth, moon, sun)
    f.step_days = 5.0

    t, y = find_minima(start_time, end_time, f, num=4)

    jd, fr = t.whole, t.tdb_fraction
    b, velocity = earth_barycenter.compute_and_differentiate(jd, fr)
    e = earth.compute(jd, fr)
    m = moon.compute(jd, fr)
    s = sun.compute(jd, fr)

    earth_to_sun = s - b - e
    moon_to_earth = e - m
    moon_to_sun = s - b - m

    earth_to_sun_AU = earth_to_sun / AU_KM
    add_aberration(earth_to_sun_AU, velocity / AU_KM,
                   length_of(earth_to_sun_AU) / C_AUDAY)
    earth_to_sun = earth_to_sun_AU * AU_KM

    solar_radius_km = 696340.0
    moon_radius_km = 1737.1
    earth_radius_km = ERAD / 1000

    angle_penumbra = arctan2(moon_radius_km + solar_radius_km, length_of(moon_to_sun))
    origin_pen_from_sun = (length_of(moon_to_sun) * solar_radius_km /
                           (moon_radius_km + solar_radius_km))
    origin_pen = s - origin_pen_from_sun * moon_to_sun / length_of(moon_to_sun)
    pen_to_earth = b + e - origin_pen
    angle_earth_from_pen = angle_between(-moon_to_sun, pen_to_earth)

    earth_from_pen_cone = (sin(angle_earth_from_pen - angle_penumbra) *
                           length_of(pen_to_earth))
    is_eclipse = ((angle_earth_from_pen < angle_penumbra) +
                  (earth_from_pen_cone <= earth_radius_km))

    origin_umbra_from_sun = (solar_radius_km * length_of(moon_to_sun) /
                             (solar_radius_km - moon_radius_km))
    angle_umbra = arctan2(solar_radius_km, origin_umbra_from_sun)
    origin_umbra = s - origin_umbra_from_sun * moon_to_sun / length_of(moon_to_sun)
    umbra_to_earth = b + e - origin_umbra

    angle_earth_from_umbra = angle_between(moon_to_sun, umbra_to_earth)

    is_angle_inside_umbra = angle_umbra - angle_earth_from_umbra > 0
    is_umbra_inside_earth_radius = length_of(umbra_to_earth) <= earth_radius_km

    angle_antumbra = pi - angle_umbra
    earth_from_antumbra_cone = sin(angle_antumbra - angle_earth_from_umbra) * length_of(
        umbra_to_earth
    )
    is_eclipse_annular = (angle_earth_from_umbra - angle_antumbra > 0) + (
            earth_from_antumbra_cone <= earth_radius_km
    )
    is_eclipse_annular *= [not x for x in is_umbra_inside_earth_radius]

    earth_from_umbra_cone = sin(angle_earth_from_umbra - angle_umbra) * length_of(
        umbra_to_earth
    )
    is_earth_radius_inside_umbra = earth_radius_km >= earth_from_umbra_cone

    is_eclipse_total = is_eclipse * (
            is_angle_inside_umbra
            + ((is_earth_radius_inside_umbra) *
               (angle_earth_from_umbra - angle_umbra < pi / 2))
            + is_umbra_inside_earth_radius
    )

    times = t[is_eclipse]
    code = is_eclipse_total.astype(byte)
    code += 2 * (is_eclipse_annular.astype(byte) -
                 (is_eclipse_annular * is_eclipse_total).astype(byte))
    code = code[is_eclipse]

    mag = empty(sum(is_eclipse))
    obs = empty(sum(is_eclipse))

    for i in range(len(times)):
        mag[i], obs[i] = _magnitude_obscurity(
            array([array(x[is_eclipse][i]) for x in earth_to_sun]),
            -array([array(x[is_eclipse][i]) for x in moon_to_earth]),
            solar_radius_km,
            moon_radius_km,
            earth_radius_km,
        )

    return times, code, mag, obs
