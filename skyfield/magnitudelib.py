# -*- coding: utf-8 -*-
"""Routines for computing magnitudes.

Planetary routines adapted from:

https://arxiv.org/pdf/1808.01973.pdf

Which links to:

https://sourceforge.net/projects/planetary-magnitudes/

Which has directories with three successive versions of their magnitude
computation, the most recent of which provides the files on which this
Python code is based:

Ap_Mag_V3.f90
Ap_Mag_Output_V3.txt
Ap_Mag_Input_V3.txt

* ``r`` planet’s distance from the Sun.
* ``delta`` from Earth?
* ``ph_ang`` illumination phase angle (degrees)

"""
from numpy import log10
from .constants import RAD2DEG
from .functions import angle_between, length_of
from .naifcodes import _target_name

def planetary_magnitude(position):
    """Given the position of a planet, return its visual magnitude.

    This prototype function — which so far only supports Mercury, Venus,
    Earth, Jupiter, and Uranus — computes the visual magnitude of a
    planet, given its position relative to an observer.

    >>> from skyfield.api import load
    >>> from skyfield.magnitudelib import planetary_magnitude
    >>> ts = load.timescale()
    >>> t = ts.utc(2020, 7, 31)
    >>> eph = load('de421.bsp')
    >>> astrometric = eph['earth'].at(t).observe(eph['jupiter barycenter'])
    >>> print('%.2f' % planetary_magnitude(astrometric))
    -2.73

    The routine does not yet take into account whether the observer is
    facing the equator or poles of Uranus, so will only be accurate to
    within about 0.1 magnitudes.

    """
    target = position.target
    function = _FUNCTIONS.get(target)
    if function is None:
        name = _target_name(target)
        raise ValueError('cannot compute the magnitude of target %s' % name)

    # Shamelessly treat the Sun as sitting at the Solar System Barycenter.
    sun_to_observer = position.center_barycentric.position.au
    observer_to_planet = position.position.au
    sun_to_planet = sun_to_observer + observer_to_planet

    r = length_of(sun_to_planet)
    delta = length_of(observer_to_planet)
    ph_ang = angle_between(-sun_to_planet, -observer_to_planet) * RAD2DEG

    return function(r, delta, ph_ang)

def _mercury_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    ph_ang_factor = (
        6.3280e-02 * ph_ang
        - 1.6336e-03 * ph_ang**2
        + 3.3644e-05 * ph_ang**3
        - 3.4265e-07 * ph_ang**4
        + 1.6893e-09 * ph_ang**5
        - 3.0334e-12 * ph_ang**6
    )
    return -0.613 + distance_mag_factor + ph_ang_factor

def _venus_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    if ph_ang < 163.7:
        ph_ang_factor = (
            -1.044E-03 * ph_ang
            + 3.687E-04 * ph_ang**2
            - 2.814E-06 * ph_ang**3
            + 8.938E-09 * ph_ang**4
        )
    else:
        ph_ang_factor = (
            236.05828 + 4.384
            - 2.81914E+00 * ph_ang
            + 8.39034E-03 * ph_ang**2
        )
    return -4.384 + distance_mag_factor + ph_ang_factor

def _earth_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10 (r * delta)
    ph_ang_factor = -1.060e-03 * ph_ang + 2.054e-04 * ph_ang**2
    return -3.99 + distance_mag_factor + ph_ang_factor

def _jupiter_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    geocentric_phase_angle_limit = 12.0

    if ph_ang <= geocentric_phase_angle_limit:  # TODO: time arrays
        ph_ang_factor = -3.7E-04 * ph_ang + 6.16E-04 * ph_ang**2
    else:
        ph_ang_factor = -2.5 * log10(
            1.0 - 1.507 * (ph_ang / 180.)
            - 0.363 * (ph_ang / 180.)**2
            - 0.062 * (ph_ang / 180.)**3
            + 2.809 * (ph_ang / 180.)**4
            - 1.876 * (ph_ang / 180.)**5
        )

    if ph_ang <= geocentric_phase_angle_limit:
        ap_mag = -9.395 + distance_mag_factor + ph_ang_factor
    else:
        ap_mag = -9.428 + distance_mag_factor + ph_ang_factor

    return ap_mag

def _uranus_magnitude(r, delta, ph_ang,
                      sun_sub_lat_planetog=0.0, earth_sub_lat_planetog=0.0):
    distance_mag_factor = 5.0 * log10 (r * delta)
    sub_lat_planetog = (abs(sun_sub_lat_planetog)
                        + abs(earth_sub_lat_planetog)) / 2.0
    sub_lat_factor = -0.00084 * sub_lat_planetog
    geocentric_phase_angle_limit = 3.1
    ap_mag = -7.110 + distance_mag_factor + sub_lat_factor
    if ph_ang > geocentric_phase_angle_limit:  # TODO: time arrays
        ap_mag += 6.587e-3 * ph_ang + 1.045e-4 * ph_ang**2
    return ap_mag

_FUNCTIONS = {
    199: _mercury_magnitude,
    299: _venus_magnitude,
    399: _earth_magnitude,
    599: _jupiter_magnitude,
    799: _uranus_magnitude,

    # Some planets can be reasonably identified with their barycenter.
    1: _mercury_magnitude,
    2: _venus_magnitude,
    5: _jupiter_magnitude,
    7: _uranus_magnitude,
}
