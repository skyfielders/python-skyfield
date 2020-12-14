from numpy import abs, einsum, sqrt, where

from .constants import C, AU_M, C_AUDAY, GS
from .functions import _AVOID_DIVIDE_BY_ZERO, dots, length_of

deflectors = ['sun', 'jupiter', 'saturn', 'moon', 'venus', 'uranus', 'neptune']
rmasses = {
    # earth-moon barycenter: 328900.561400
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

    Given the ICRS `position` [x,y,z] of an object (au) that is being
    viewed from the `observer` also expressed as [x,y,z], and given an
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

    jd_tdb = t.tdb
    ts = t.ts
    for name in deflectors[:count]:
        try:
            deflector = ephemeris[name]
        except KeyError:
            deflector = ephemeris[name + ' barycenter']

        # Get position of gravitating body wrt ss barycenter at time 't_tdb'.

        bposition = deflector.at(ts.tdb(jd=jd_tdb)).position.au  # TODO

        # Get position of gravitating body wrt observer at time 'jd_tdb'.

        gpv = bposition - observer

        # Compute light-time from point on incoming light ray that is closest
        # to gravitating body.

        dlt = light_time_difference(position, gpv)

        # Get position of gravitating body wrt ss barycenter at time when
        # incoming photons were closest to it.

        tclose = jd_tdb

        # if dlt > 0.0:
        #     tclose = jd - dlt

        tclose = where(dlt > 0.0, jd_tdb - dlt, tclose)
        tclose = where(tlt < dlt, jd_tdb - tlt, tclose)

        # if tlt < dlt:
        #     tclose = jd - tlt

        bposition = deflector.at(ts.tdb(jd=tclose)).position.au  # TODO
        rmass = rmasses[name]
        _add_deflection(position, observer, bposition, rmass)

    # If observer is not at geocenter, add in deflection due to Earth.

    if include_earth_deflection.any():
        deflector = ephemeris['earth']
        bposition = deflector.at(ts.tdb(jd=tclose)).position.au  # TODO
        rmass = rmasses['earth']
        # TODO: Make the following code less messy, maybe by having
        # _add_deflection() return a new vector instead of modifying the
        # old one in-place.
        deflected_position = position.copy()
        _add_deflection(deflected_position, observer, bposition, rmass)
        if include_earth_deflection.shape:
            position[:,include_earth_deflection] = (
                deflected_position[:,include_earth_deflection])
        else:
            position[:] = deflected_position[:]

def light_time_difference(position, observer_position):
    """Returns the difference in light-time, for a star,
      between the barycenter of the solar system and the observer (or
      the geocenter).

    """
    # From 'pos1', form unit vector 'u1' in direction of star or light
    # source.

    dis = length_of(position)
    u1 = position / (dis + _AVOID_DIVIDE_BY_ZERO)

    # Light-time returned is the projection of vector 'pos_obs' onto the
    # unit vector 'u1' (formed from 'pos1'), divided by the speed of light.

    diflt = einsum('a...,a...', u1, observer_position) / C_AUDAY
    return diflt

def _add_deflection(position, observer, deflector, rmass):
    """Correct a position vector for how one particular mass deflects light.

    Given the ICRS `position` [x,y,z] of an object (AU) together with
    the positions of an `observer` and a `deflector` of reciprocal mass
    `rmass`, this function updates `position` in-place to show how much
    the presence of the deflector will deflect the image of the object.

    """
    # Construct vector 'pq' from gravitating body to observed object and
    # construct vector 'pe' from gravitating body to observer.

    pq = observer + position - deflector
    pe = observer - deflector

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
    # is set to zero set 'pos2' equal to 'pos1'.

    make_no_correction = abs(edotp) > 0.99999999999

    # Compute scalar factors.

    fac1 = 2.0 * GS / (C * C * emag * AU_M * rmass)
    fac2 = 1.0 + qdote

    # Correct position vector.

    position += where(make_no_correction, 0.0,
                      fac1 * (pdotq * ehat - edotp * qhat) / fac2 * pmag)

def add_aberration(position, velocity, light_time):
    """Correct a relative position vector for aberration of light.

    Given the relative `position` [x,y,z] of an object (AU) from a
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
