from numpy import abs, clip, einsum, sqrt, where

from .constants import C, AU_M, C_AUDAY, GS
from .functions import _AVOID_DIVIDE_BY_ZERO, dots, length_of

deflectors = ['sun', 'jupiter', 'saturn', 'moon', 'venus', 'uranus', 'neptune']
rmasses = {
    # earth-moon barycenter: 328900.561400
    199: 6023600.0,             # mercury
    299: 408523.71,             # venus
    399: 332946.050895,         # earth
    499: 3098708.0,             # mars
    599: 1047.3486,             # jupiter
    699: 3497.898,              # saturn
    799: 22902.98,              # uranus
    899: 19412.24,              # neptune
    999: 135200000.0,           # pluto
    10: 1.0,                    # sun
    301: 27068700.387534,       # moon

    # For compatibility with any user code that discovered and used this
    # undocumented `rmasses` dict prior to Skyfield 1.53:
    'mercury': 6023600.0,
    'venus': 408523.71,
    'earth': 332946.050895,
    'mars': 3098708.0,
    'jupiter': 1047.3486,
    'saturn': 3497.898,
    'uranus': 22902.98,
    'neptune': 19412.24,
    'pluto': 135200000.0,
    'sun': 1.0,
    'moon': 27068700.387534,
}

def add_deflection(position, observer, ephemeris, t,
                   include_earth_deflection, count=3):
    """Update `position` for how solar system masses will deflect its light.

    Given the ICRS `position` |xyz| of an object (au) that is being
    viewed from the `observer` also expressed as |xyz|, and given an
    ephemeris that can be used to determine solar system body positions,
    and given the time `t` and Boolean `apply_earth` indicating whether
    to worry about the effect of Earth's mass, and a `count` of how many
    major solar system bodies to worry about, this function updates
    `position` in-place to show how the masses in the solar system will
    deflect its image.

    """
    # Compute light-time to observed object.

    tlt = length_of(position) / C_AUDAY

    # Cycle through gravitating bodies.

    for name in deflectors[:count]:
        try:
            deflector = ephemeris[name]
        except KeyError:
            deflector = ephemeris[name + ' barycenter']

        rmass = rmasses[name]
        pe = _compute_deflector_position(
            t, observer, position, deflector, tlt,
        )
        position += _compute_deflection(position, pe, rmass)

    # If observer is not at geocenter, add in deflection due to Earth.

    if include_earth_deflection.any():
        deflector = ephemeris['earth']
        bposition = deflector.at(t).xyz.au
        rmass = rmasses['earth']
        pe = observer - bposition
        d = _compute_deflection(position, pe, rmass)
        if include_earth_deflection.shape:
            d *= include_earth_deflection  # where False, set `d` to zero
        position += d

def _compute_deflector_position(t, observer, position, deflector, tlt):
    """Compute where a deflector was when closest to a light beam."""

    # Get position of gravitating body wrt observer.

    bposition = deflector.at(t).xyz.au
    gpv = bposition - observer

    # Compute light-time from point on incoming light ray that is closest
    # to gravitating body.

    dlt = light_time_difference(position, gpv)

    # Should really add TDB offset, not TT, but doesn't matter.

    tclose = t - clip(dlt, 0.0, tlt)

    # Get position of gravitating body wrt ss barycenter at time when
    # incoming photons were closest to it.

    bposition = deflector.at(tclose).xyz.au
    pe = observer - bposition
    return pe

def light_time_difference(position, deflector_position):
    """When did the light from ``position`` pass closest to a given deflector?

    Given an observer at the origin and a ``position`` |xyz| in AU, how
    recently did the light from the position pass closest to the object
    whose |xyz| is given as ``deflector_position``?  The answer is
    returned as a floating point number of days.

    """
    # From 'pos1', form unit vector 'u1' in direction of star or light
    # source.

    dis = length_of(position)
    u1 = position / (dis + _AVOID_DIVIDE_BY_ZERO)

    # Light-time returned is the projection of vector 'pos_obs' onto the
    # unit vector 'u1' (formed from 'pos1'), divided by the speed of light.

    diflt = einsum('a...,a...', u1, deflector_position) / C_AUDAY
    return diflt

def _compute_deflection(position, pe, rmass):
    """Compute how much a mass will deflect a position.

    The ``position`` is relative to the observer in AU; ``pe`` is the
    position of the observer relative to the deflector in AU; and
    ``rmass`` is the deflector's reciprocal mass.

    """
    pq = position + pe

    # Compute vector magnitudes and unit vectors.

    pmag = length_of(position)
    qmag = length_of(pq)
    emag = length_of(pe)

    phat = position / where(pmag, pmag, 1.0)  # where() avoids divide-by-zero
    qhat = pq / where(qmag, qmag, 1.0)
    ehat = pe / where(emag, emag, 1.0)

    # Compute dot products of vectors.

    pdotq = dots(phat, qhat)
    qdote = dots(qhat, ehat)
    edotp = dots(ehat, phat)

    # If gravitating body is observed object, or is on a straight line
    # toward or away from observed object to within 1 arcsec, deflection
    # is set to zero.

    flag = abs(edotp) <= 0.99999999999

    # Compute scalar factors.

    fac1 = 2.0 * GS / (C * C * emag * AU_M * rmass)
    fac2 = 1.0 + qdote

    return flag * fac1 * (pdotq * ehat - edotp * qhat) / fac2 * pmag

def add_aberration(position, velocity, light_time):
    """Correct a relative position vector for aberration of light.

    Given the relative `position` |xyz| of an object (AU) from a
    particular observer, the `velocity` [dx,dy,dz] at which the observer
    is traveling (AU/day), and the light propagation delay `light_time`
    to the object (days), this function updates `position` in-place to
    give the object's apparent position due to the aberration of light.

    """
    p1mag = light_time * C_AUDAY
    vemag = length_of(velocity)
    beta = vemag / C_AUDAY
    dot = dots(position, velocity)

    cosd = dot / (p1mag * vemag + _AVOID_DIVIDE_BY_ZERO)
    gammai = sqrt(1.0 - beta * beta)
    p = beta * cosd
    q = (1.0 + p / (1.0 + gammai)) * light_time
    r = 1.0 + p

    position *= gammai
    position += q * velocity
    position /= r
