"""Formulae for specific earth behaviors and effects."""

from numpy import (arcsin, arccos, array, clip, cos, einsum, fmod,
                   minimum, pi, sin, sqrt, zeros_like)

from .constants import (AU, ANGVEL, AU_KM, DAY_S, DEG2RAD, ERAD,
                        ERAD_KM, IERS_2010_INVERSE_EARTH_FLATTENING,
                        RAD2DEG, T0)
from .framelib import J2000_to_ICRS
from .functions import dots
from .nutationlib import earth_tilt

rade = ERAD / AU
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

    pos, vel = terra(topos, gast)

    t = einsum('ji...,kj...->ik...', jd.N, jd.P)
    t = einsum('ij...,jk...->ik...', t, J2000_to_ICRS)

    pos = einsum('i...,ij...->j...', pos, t)
    vel = einsum('i...,ij...->j...', vel, t)

    return pos, vel


def terra(topos, st):
    """Compute the position and velocity of a terrestrial observer.

    `topos` - `Topos` object describing a geographic position.
    `st` - Array of sidereal times in floating-point hours.

    The return value is a tuple of two 3-vectors `(pos, vel)` in the
    dynamical reference system whose components are measured in AU with
    respect to the center of the Earth.

    """
    zero = zeros_like(st)
    phi = topos.latitude
    sinphi = sin(phi)
    cosphi = cos(phi)
    c = 1.0 / sqrt(cosphi * cosphi +
                   sinphi * sinphi * one_minus_flattening_squared)
    s = one_minus_flattening_squared * c
    ht_km = topos.elevation / 1000.0
    ach = ERAD_KM * c + ht_km
    ash = ERAD_KM * s + ht_km

    # Compute local sidereal time factors at the observer's longitude.

    stlocl = st * 15.0 * DEG2RAD + topos.longitude
    sinst = sin(stlocl)
    cosst = cos(stlocl)

    # Compute position vector components in kilometers.

    ac = ach * cosphi
    pos = array((ac * cosst, ac * sinst, zero + ash * sinphi)) / AU_KM

    # Compute velocity vector components in kilometers/sec.

    aac = ANGVEL * ach * cosphi
    vel = array((-aac * sinst, aac * cosst, zero)) / AU_KM * DAY_S

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
