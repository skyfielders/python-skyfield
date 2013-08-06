from numpy import abs, sqrt, where
from skyfield.functions import dots, length

C = 299792458.0
AU = 1.4959787069098932e+11
C_AUDAY = 173.1446326846693

# Heliocentric gravitational constant in meters^3 / second^2, from DE-405.

GS = 1.32712440017987e+20

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

def add_deflection(position, observer, ephemeris, jd_tdb,
                   include_earth_deflection, count=3):
    """Update `position` for how solar system masses will deflect its light.

    Given the ICRS `position` [x,y,z] of an object (AU) that is being
    viewed from the `observer` also expressed as [x,y,z], and given an
    ephemeris that can be used to determine solar system body positions,
    and given the time `jd` and Boolean `apply_earth` indicating whether
    to worry about the effect of Earth's mass, and a `count` of how many
    major solar system bodies to worry about, this function updates
    `position` in-place to show how the masses in the solar system will
    deflect its image.

    """
    # Compute light-time to observed object.

    tlt = length(position) / C_AUDAY

    # Cycle through gravitating bodies.

    for name in deflectors[:count]:

        # Get position of gravitating body wrt ss barycenter at time 'jd_tdb'.

        bpv = ephemeris.compute(name, jd_tdb)

        # Get position of gravitating body wrt observer at time 'jd_tdb'.

        gpv = bpv.position - observer

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

        bpv = ephemeris.compute(name, tclose)
        rmass = rmasses[name]
        _add_deflection(position, observer, bpv.position, rmass)

    # If observer is not at geocenter, add in deflection due to Earth.

    if include_earth_deflection.any():
        bpv = ephemeris.compute('earth', jd_tdb)
        rmass = rmasses['earth']
        _add_deflection(position, observer, bpv.position, rmass)

#

def light_time_difference(position, observer_position):
    """Returns the difference in light-time, for a star,
      between the barycenter of the solar system and the observer (or
      the geocenter).

    """
    # From 'pos1', form unit vector 'u1' in direction of star or light
    # source.

    dis = length(position)
    u1 = position / dis

    # Light-time returned is the projection of vector 'pos_obs' onto the
    # unit vector 'u1' (formed from 'pos1'), divided by the speed of light.

    diflt = (observer_position * u1).sum(axis=0) / C_AUDAY
    return diflt

#

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

    pmag = length(position)
    qmag = length(pq)
    emag = length(pe)

    phat = position / pmag
    qhat = pq / qmag
    ehat = pe / emag

    # Compute dot products of vectors.

    pdotq = dots(phat, qhat)
    qdote = dots(qhat, ehat)
    edotp = dots(ehat, phat)

    # If gravitating body is observed object, or is on a straight line
    # toward or away from observed object to within 1 arcsec, deflection
    # is set to zero set 'pos2' equal to 'pos1'.

    make_no_correction = abs(edotp) > 0.99999999999

    # Compute scalar factors.

    fac1 = 2.0 * GS / (C * C * emag * AU * rmass)
    fac2 = 1.0 + qdote

    # Correct position vector.

    position += where(make_no_correction, 0.0,
                      fac1 * (pdotq * ehat - edotp * qhat) / fac2 * pmag)

#

def add_aberration(position, velocity, lighttime):
    """Correct a relative position vector for aberration of light.

    Given the relative `position` [x,y,z] of an object (AU) from a
    particular observer, the `velocity` [dx,dy,dz] at which the observer
    is traveling (AU/day), and the light propagation delay `lighttime`
    to the object (days), this function updates `position` in-place to
    give the object's apparent position due to the aberration of light.

    """
    p1mag = lighttime * C_AUDAY
    vemag = length(velocity)
    beta = vemag / C_AUDAY
    dot = dots(position, velocity)

    cosd = dot / (p1mag * vemag)
    gammai = sqrt(1.0 - beta * beta)
    p = beta * cosd
    q = (1.0 + p / (1.0 + gammai)) * lighttime
    r = 1.0 + p

    position *= gammai
    position += q * velocity
    position /= r
