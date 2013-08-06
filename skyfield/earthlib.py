"""Formulae for specific earth behaviors and effects."""

from numpy import (arcsin, arccos, array, clip, cos, einsum,
                   pi, sin, sqrt, zeros_like)
from skyfield import timescales
from skyfield.angles import DEG2RAD
from skyfield.framelib import J2000_to_ICRS
from skyfield.functions import dots
from skyfield.nutationlib import earth_tilt, compute_nutation
from skyfield.precessionlib import compute_precession

ANGVEL = 7.2921150e-5
AU = 1.4959787069098932e+11
AU_KM = 1.4959787069098932e+8
ERAD = 6378136.6
ERAD_KM = ERAD / 1000.0
F = 0.003352819697896
RAD2DEG = 57.295779513082321
rade = ERAD / AU
halfpi = pi / 2.0


def geocentric_position_and_velocity(topos, jd_tt, delta_t=0.0):
    """Compute the geocentric position, velocity of a terrestrial observer.

    `topos` - `Topos` object describing a location.
    `jd_tt` - Array of Julian dates in Terrestrial Time.

    The return value is a 2-element tuple `(pos, vel)` of 3-vectors
    which each measure position in AU long the axes of the ICRS.

    """
    jd_tdb = jd_tt + timescales.tdb_minus_tt(jd_tt)
    jd_ut1 = jd_tt - (delta_t / 86400.)

    gmst = timescales.sidereal_time(jd_ut1, delta_t)
    x1, x2, eqeq, x3, x4 = earth_tilt(jd_tdb)
    gast = gmst + eqeq / 3600.0

    pos, vel = terra(topos, gast)

    n = compute_nutation(jd_tdb)
    p = compute_precession(jd_tdb)
    f = J2000_to_ICRS

    t = einsum('jin,kjn->ikn', n, p)
    t = einsum('ijn,jk->ikn', t, f)

    pos = einsum('in,ijn->jn', pos, t)
    vel = einsum('in,ijn->jn', vel, t)

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
    df = 1.0 - F
    df2 = df * df

    phi = topos.latitude
    sinphi = sin(phi)
    cosphi = cos(phi)
    c = 1.0 / sqrt(cosphi * cosphi + df2 * sinphi * sinphi)
    s = df2 * c
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
    vel = array((-aac * sinst, aac * cosst, zero)) / AU_KM * 86400.0

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

    if disobs >= rade:
        aprad = arcsin(rade / disobs)
    else:
        aprad = halfpi

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
