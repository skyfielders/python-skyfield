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
from numpy import log10, where
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
    condition = ph_ang < 163.7
    a0 = where(condition, 0.0, 236.05828 + 4.384)
    a1 = where(condition, - 1.044E-03, - 2.81914E+00)
    a2 = where(condition, + 3.687E-04, + 8.39034E-03)
    a3 = where(condition, - 2.814E-06, 0.0)
    a4 = where(condition, + 8.938E-09, 0.0)
    ph_ang_factor = a4
    for a in a3, a2, a1, a0:
        ph_ang_factor *= ph_ang
        ph_ang_factor += a
    return -4.384 + distance_mag_factor + ph_ang_factor

def _earth_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    ph_ang_factor = -1.060e-03 * ph_ang + 2.054e-04 * ph_ang**2
    return -3.99 + distance_mag_factor + ph_ang_factor

def _mars_magnitude(r, delta, ph_ang):
    r_mag_factor = 2.5 * log10(r * r)
    delta_mag_factor = 2.5 * log10(delta * delta)
    distance_mag_factor = r_mag_factor + delta_mag_factor

    geocentric_phase_angle_limit = 50.0

    condition = ph_ang <= geocentric_phase_angle_limit
    a = where(condition, 2.267E-02, - 0.02573)
    b = where(condition, - 1.302E-04, 0.0003445)
    ph_ang_factor = a * ph_ang + b * ph_ang**2

    # Compute the effective central meridian longitude
    # eff_CM = ( sub_earth_long + sub_sun_long ) / 2.
    # if ( abs ( sub_earth_long - sub_sun_long ) > 180. ):
    #     Eff_CM = Eff_CM + 180.
    # if ( Eff_CM > 360. ):
    #     Eff_CM = Eff_CM - 360.

    # ! Use Stirling interpolation to determine the magnitude correction
    # call Mars_Stirling ( 'R', eff_CM, mag_corr_rot )

    # Convert the ecliptic longitude to Ls
    # Ls = h_ecl_long + Ls_offset
    # if ( Ls > 360. ) Ls = Ls - 360.
    # if ( Ls <   0. ) Ls = Ls + 360.

    # Use Stirling interpolation to determine the magnitude correction
    # call Mars_Stirling ( 'O', Ls, mag_corr_orb )

    # Until effects from Mars rotation are written up:
    mag_corr_rot = 0.0
    mag_corr_orb = 0.0

    # Add factors to determine the apparent magnitude
    ap_mag = where(ph_ang <= geocentric_phase_angle_limit, -1.601, -0.367)
    ap_mag += distance_mag_factor + ph_ang_factor + mag_corr_rot + mag_corr_orb

    return ap_mag

def _jupiter_magnitude(r, delta, ph_ang):
    distance_mag_factor = 5 * log10(r * delta)
    geocentric_phase_angle_limit = 12.0
    ph_ang_pi = ph_ang / 180.0
    ph_ang_factor = where(
        ph_ang <= geocentric_phase_angle_limit,
        (6.16E-04 * ph_ang - 3.7E-04) * ph_ang,
        -2.5 * log10(
            ((((- 1.876 * ph_ang_pi + 2.809) * ph_ang_pi - 0.062) * ph_ang_pi
              - 0.363) * ph_ang_pi - 1.507) * ph_ang_pi + 1.0
        ),
    )
    ap_mag = where(
        ph_ang <= geocentric_phase_angle_limit,
        -9.395 + distance_mag_factor + ph_ang_factor,
        -9.428 + distance_mag_factor + ph_ang_factor,
    )
    return ap_mag

def _uranus_magnitude(r, delta, ph_ang,
                      sun_sub_lat_planetog=0.0, earth_sub_lat_planetog=0.0):
    distance_mag_factor = 5.0 * log10 (r * delta)
    sub_lat_planetog = (abs(sun_sub_lat_planetog)
                        + abs(earth_sub_lat_planetog)) / 2.0
    sub_lat_factor = -0.00084 * sub_lat_planetog
    geocentric_phase_angle_limit = 3.1
    ap_mag = -7.110 + distance_mag_factor + sub_lat_factor
    ap_mag += where(
        ph_ang > geocentric_phase_angle_limit,
        (1.045e-4 * ph_ang + 6.587e-3) * ph_ang,
        0.0,
    )
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
