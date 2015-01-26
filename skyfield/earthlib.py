"""Formulae for specific earth behaviors and effects."""

from numpy import (abs, arcsin, arccos, array, clip, cos, einsum, exp, fmod,
                   logical_and, minimum, pi, sin, sqrt, tan, where, zeros_like)

from .constants import (AU_M, ANGVEL, DAY_S, DEG2RAD, ERAD,
                        IERS_2010_INVERSE_EARTH_FLATTENING, RAD2DEG, T0)
from .functions import dots
from .nutationlib import earth_tilt

rade = ERAD / AU_M
one_minus_flattening = 1.0 - 1.0 / IERS_2010_INVERSE_EARTH_FLATTENING
one_minus_flattening_squared = one_minus_flattening * one_minus_flattening


def geocentric_position_and_velocity(topos, jd):
    """Compute the GCRS position and velocity of a terrestrial observer.

    `topos` - `Topos` object describing a location.
    `jd` - a JulianDate.

    The return value is a 2-element tuple `(pos, vel)` of 3-vectors
    which each measure position in AU long the axes of the ICRS.

    """
    gmst = sidereal_time(jd)
    x1, x2, eqeq, x3, x4 = earth_tilt(jd)
    gast = gmst + eqeq / 3600.0

    pos, vel = terra(
        topos.latitude.radians,
        topos.longitude.radians,
        topos.elevation.m,
        gast,
        )

    pos = einsum('ij...,j...->i...', jd.MT, pos)
    vel = einsum('ij...,j...->i...', jd.MT, vel)

    return pos, vel


def terra(latitude, longitude, elevation, st):
    """Compute the position and velocity of a terrestrial observer.

    `latitude` - Latitude in radians.
    `longitude` - Longitude in radians.
    `elevation` - Elevation in meters.
    `st` - Array of sidereal times in floating-point hours.

    The return value is a tuple of two 3-vectors `(pos, vel)` in the
    dynamical reference system whose components are measured in AU with
    respect to the center of the Earth.

    """
    zero = zeros_like(st)
    sinphi = sin(latitude)
    cosphi = cos(latitude)
    c = 1.0 / sqrt(cosphi * cosphi +
                   sinphi * sinphi * one_minus_flattening_squared)
    s = one_minus_flattening_squared * c
    ht = elevation
    ach = ERAD * c + ht
    ash = ERAD * s + ht

    # Compute local sidereal time factors at the observer's longitude.

    stlocl = st * 15.0 * DEG2RAD + longitude
    sinst = sin(stlocl)
    cosst = cos(stlocl)

    # Compute position vector components in kilometers.

    ac = ach * cosphi
    pos = array((ac * cosst, ac * sinst, zero + ash * sinphi)) / AU_M

    # Compute velocity vector components in kilometers/sec.

    aac = ANGVEL * ach * cosphi
    vel = array((-aac * sinst, aac * cosst, zero)) / AU_M * DAY_S

    return pos, vel


def compute_limb_angle(position, observer):
    """Determine the angle of an object above or below the Earth's limb.

    Given an object's GCRS `position` [x,y,z] in AU and the position of
    an `observer` in the same coordinate system, return a tuple that is
    composed of `(limb_ang, nadir_ang)`:

    limb_angle
        Angle of observed object above (+) or below (-) limb in degrees.
    nadir_angle
        Nadir angle of observed object as a fraction of apparent radius
        of limb: <1.0 means below the limb, =1.0 means on the limb, and
        >1.0 means above the limb.

    """
    # Compute the distance to the object and the distance to the observer.

    disobj = sqrt(dots(position, position))
    disobs = sqrt(dots(observer, observer))

    # Compute apparent angular radius of Earth's limb.

    aprad = arcsin(minimum(rade / disobs, 1.0))

    # Compute zenith distance of Earth's limb.

    zdlim = pi - aprad

    # Compute zenith distance of observed object.

    coszd = dots(position, observer) / (disobj * disobs)
    coszd = clip(coszd, -1.0, 1.0)
    zdobj = arccos(coszd)

    # Angle of object wrt limb is difference in zenith distances.

    limb_angle = (zdlim - zdobj) * RAD2DEG

    # Nadir angle of object as a fraction of angular radius of limb.

    nadir_angle = (pi - zdobj) / aprad

    return limb_angle, nadir_angle


def sidereal_time(jd, use_eqeq=False):
    """Compute Greenwich sidereal time at Julian date `jd_ut1`."""

    t = (jd.tdb - T0) / 36525.0

    # Equation of equinoxes.

    if use_eqeq:
        ee = earth_tilt(jd)[2]
        eqeq = ee * 15.0
    else:
        eqeq = 0.0

    # Compute the Earth Rotation Angle.  Time argument is UT1.

    theta = earth_rotation_angle(jd.ut1)

    # The equinox method.  See Circular 179, Section 2.6.2.
    # Precession-in-RA terms in mean sidereal time taken from third
    # reference, eq. (42), with coefficients in arcseconds.

    st = eqeq + ( 0.014506 +
        (((( -    0.0000000368   * t
             -    0.000029956  ) * t
             -    0.00000044   ) * t
             +    1.3915817    ) * t
             + 4612.156534     ) * t)

    # Form the Greenwich sidereal time.

    gst = fmod((st / 3600.0 + theta), 360.0) / 15.0

    gst += 24.0 * (gst < 0.0)

    return gst


def earth_rotation_angle(jd_ut1):
    """Return the value of the Earth Rotation Angle (theta) for a UT1 date.

    Uses the expression from the note to IAU Resolution B1.8 of 2000.

    """
    thet1 = 0.7790572732640 + 0.00273781191135448 * (jd_ut1 - T0)
    thet3 = jd_ut1 % 1.0
    return (thet1 + thet3) % 1.0 * 360.0


def refraction(alt_degrees, temperature_C, pressure_mbar):
    """Given an observed altitude, return how much the image is refracted.

    Zero refraction is returned both for objects very near the zenith,
    as well as for objects more than one degree below the horizon.

    """
    r = 0.016667 / tan((alt_degrees + 7.31 / (alt_degrees + 4.4)) * DEG2RAD)
    d = r * (0.28 * pressure_mbar / (temperature_C + 273.0))
    return where((-1.0 <= alt_degrees) & (alt_degrees <= 89.9), d, 0.0)


def refract(alt_degrees, temperature_C, pressure_mbar):
    """Given an unrefracted `alt` determine where it will appear in the sky."""
    alt = alt_degrees
    while True:
        alt1 = alt
        alt = alt_degrees + refraction(alt, temperature_C, pressure_mbar)
        converged = abs(alt - alt1) <= 3.0e-5
        if converged.all():
            break
    return alt
