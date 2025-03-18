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
from numpy import array, clip, exp, log10, nan, sin, where
from .constants import RAD2DEG
from .functions import angle_between, length_of
from .naifcodes import _target_name

# See "design/planet_tilts.py" in the Skyfield repository.
_SATURN_POLE = array([0.08547883, 0.07323576, 0.99364475])
_SATURN_POLE_2D = _SATURN_POLE[:, None]
_URANUS_POLE = array([-0.21199958, -0.94155916, -0.26176809])
_URANUS_POLE_2D = _URANUS_POLE[:, None]

def planetary_magnitude(position):
    """Given the position of a planet, return its visual magnitude.

    >>> from skyfield.api import load
    >>> from skyfield.magnitudelib import planetary_magnitude
    >>> ts = load.timescale()
    >>> t = ts.utc(2020, 7, 31)
    >>> eph = load('de421.bsp')
    >>> astrometric = eph['earth'].at(t).observe(eph['jupiter barycenter'])
    >>> print('%.2f' % planetary_magnitude(astrometric))
    -2.73

    The formulae are from `Mallama and Hilton “Computing Apparent
    Planetary Magnitude for the Astronomical Almanac” (2018)
    <https://arxiv.org/pdf/1808.01973.pdf>`_.  Two of the formulae have
    inherent limits:

    * Saturn’s magnitude is unknown and the function will return ``nan``
      (the floating-point value “Not a Number”) if the “illumination
      phase angle” — the angle of the vertex observer-Saturn-Sun —
      exceeds 6.5°.

    * Neptune’s magnitude is unknown and will return ``nan`` if the
      illumination phase angle exceeds 1.9° and the position's date is
      before the year 2000.

    And one formula is not fully implemented (though contributions are
    welcome!):

    * Skyfield does not compute which features on Mars are facing the
      observer, which can introduce an error of ±0.06 magnitude.

    """
    target = position.target
    function = _FUNCTIONS.get(target)
    if function is None:
        name = _target_name(target)
        raise ValueError('cannot compute the magnitude of target %s' % name)

    # Shamelessly treat the Sun as sitting at the Solar System Barycenter.
    sun_to_observer = position.center_barycentric.xyz.au
    observer_to_planet = position.xyz.au
    sun_to_planet = sun_to_observer + observer_to_planet

    r = length_of(sun_to_planet)
    delta = length_of(observer_to_planet)
    ph_ang = angle_between(sun_to_planet, observer_to_planet) * RAD2DEG

    if function is _saturn_magnitude:
        if len(sun_to_planet.shape) > 1:
            pole = _SATURN_POLE_2D
        else:
            pole = _SATURN_POLE

        a = angle_between(pole, sun_to_planet)
        sun_sub_lat = a * RAD2DEG - 90.0

        a = angle_between(pole, observer_to_planet)
        observer_sub_lat = a * RAD2DEG - 90.0

        return function(r, delta, ph_ang, sun_sub_lat, observer_sub_lat)

    if function is _uranus_magnitude:
        if len(sun_to_planet.shape) > 1:
            pole = _URANUS_POLE_2D
        else:
            pole = _URANUS_POLE

        a = angle_between(pole, sun_to_planet)
        sun_sub_lat = a * RAD2DEG - 90.0

        a = angle_between(pole, observer_to_planet)
        observer_sub_lat = a * RAD2DEG - 90.0

        return function(r, delta, ph_ang, sun_sub_lat, observer_sub_lat)

    if function is _neptune_magnitude:
        year = position.t.J
        return function(r, delta, ph_ang, year)

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

def _saturn_magnitude(r, delta, ph_ang, sun_sub_lat, earth_sub_lat, rings=True):
    # Note that sun_sub_lat and earth_sub_lat should be saturnicentric
    # latitude, not saturnidetic.

    r_mag_factor = 2.5 * log10(r * r)
    delta_mag_factor = 2.5 * log10(delta * delta)
    distance_mag_factor = r_mag_factor + delta_mag_factor

    # Then take the square root of the product of the saturnicentric
    # latitude of the Sun and that of Earth but set to zero when the
    # signs are opposite
    product = sun_sub_lat * earth_sub_lat
    signs_same = product >= 0.0
    square_root = product ** where(signs_same, 0.5, 0.0)  # avoid sqrt(neg)
    sub_lat_geoc = where(signs_same, square_root, 0.0)

    # Compute the effect of phase angle and inclination

    geocentric_phase_angle_limit = 6.5
    geocentric_inclination_limit = 27.0

    is_within_geocentric_bounds = (
        (ph_ang <= geocentric_phase_angle_limit)
        & (sub_lat_geoc <= geocentric_inclination_limit)
    )

    ap_mag = where(
        is_within_geocentric_bounds,
        where(
            rings,

            # Use equation #10 for globe+rings and geocentric circumstances.
            -8.914 - 1.825 * sin(sub_lat_geoc / RAD2DEG) + 0.026 * ph_ang
            - 0.378 * sin(sub_lat_geoc / RAD2DEG) * exp(-2.25 * ph_ang),

            # Use equation #11 for globe-alone and geocentric circumstances
            -8.95 - 3.7e-4 * ph_ang + 6.16e-4 * ph_ang**2,
        ),
        where(
            (ph_ang > geocentric_phase_angle_limit) & _not(rings),

            # Use equation #12 for globe-alone beyond geocentric phase
            # angle limit
            -8.94 + 2.446e-4 * ph_ang + 2.672e-4 * ph_ang**2
            - 1.506e-6 * ph_ang**3 +4.767e-9 * ph_ang**4,

            nan,
        ),
    )

    ap_mag = ap_mag + distance_mag_factor
    return ap_mag

def _not(b):
    return 1 - b  # since ~b raises a DeprecationWarning in Python 3.13

def _uranus_magnitude(r, delta, ph_ang,
                      sun_sub_lat_planetog, earth_sub_lat_planetog):
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

def _neptune_magnitude(r, delta, ph_ang, year):
    r_mag_factor = 2.5 * log10(r * r)
    delta_mag_factor = 2.5 * log10(delta * delta)
    distance_mag_factor = r_mag_factor + delta_mag_factor

    # Equation 16 compute the magnitude at unit distance as a function of time
    ap_mag = clip(-6.89 - 0.0054 * (year - 1980.0), -7.00, -6.89)
    ap_mag += distance_mag_factor

    geocentric_phase_angle_limit = 1.9

    ap_mag = where(
        ph_ang > geocentric_phase_angle_limit,

        # Add phase angle factor from equation 17
        # Check the year because equation 17 only pertains to t > 2000.0
        where(
            year >= 2000.0,
            ap_mag + 7.944e-3 * ph_ang + 9.617e-5 * ph_ang**2,
            nan,
        ),

        # Otherwise leave the value unchanged.
        ap_mag,
    )

    return ap_mag

_FUNCTIONS = {
    199: _mercury_magnitude,
    299: _venus_magnitude,
    399: _earth_magnitude,
    499: _mars_magnitude,
    599: _jupiter_magnitude,
    699: _saturn_magnitude,
    799: _uranus_magnitude,
    899: _neptune_magnitude,

    # Some planets can be reasonably identified with their barycenter.
    1: _mercury_magnitude,
    2: _venus_magnitude,
    4: _mars_magnitude,
    5: _jupiter_magnitude,
    6: _saturn_magnitude,
    7: _uranus_magnitude,
    8: _neptune_magnitude,
}
