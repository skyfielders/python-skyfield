"""An interface between Skyfield and the Python ``sgp4`` library."""

from numpy import (
    array, concatenate, cross, einsum, identity, ones_like, repeat, zeros_like
)
from sgp4.api import SGP4_ERRORS, Satrec

from .almanac import _find_discrete, _find_maxima
from .constants import AU_KM, DAY_S, T0, tau
from .functions import rot_x, rot_y, rot_z
from .positionlib import ITRF_to_GCRS2
from .timelib import Timescale
from .vectorlib import VectorFunction

_minutes_per_day = 1440.

# Since satellite calculations are done entirely in UTC, we can display
# and return their epoch using a null Timescale.  Let's just hope that
# no one ever pulls the .epoch off of an EarthSatellite and tries to
# generate high precision positions with it.  If that becomes a problem,
# we can make a more severe change to the satellite API to require a
# Timescale when the user asks for the epoch as a Time.  (It was
# probably a mistake to have an .epoch convenience attribute.)

_identity = identity(3)
_infs = array(('-inf', 'inf'), float)
_ts = Timescale(array((_infs, (0.0, 0.0))), _infs, array((37.0, 37.0)))

class EarthSatellite(VectorFunction):
    """An Earth satellite loaded from a TLE file and propagated with SGP4.

    An earth satellite object is a Skyfield vector function, so call its
    :meth:`~skyfield.vectorlib.VectorSum.at()` method to generate its
    position in the sky, or use addition and subtraction to combine it
    with other vectors.

    Satellite parameters are generally only accurate for a week or two
    around the *epoch* of the parameters, the date for which they were
    generated, which is available as an attribute:

    ``epoch``
        A Skyfield :class:`~skyfield.timelib.Time` giving the exact
        epoch moment for these satellite orbit parameters.

    When building a satellite, use the arguments ``line1`` and ``line2``
    to provide the two data lines from a TLE file as separate strings.
    Optional ``name`` lets you give a name to the satellite, accessible
    later through the ``name`` attribute.  ``ts`` is a
    :class:`~skyfield.timelib.Timescale` object, used to generate the
    ``epoch`` value; if it is not provided, the satellite will use a
    built in ``Timescale`` object.

    If you are interested in the catalog entry details, the SGP4 model
    parameters for a particular satellite can be accessed through its
    ``model`` attribute:

    ``model.satnum``
        The unique satellite NORAD catalog number given in the TLE file.
    ``model.epochyr``
        Full four-digit year of this element set's epoch moment.
    ``model.epochdays``
        Fractional days into the year of the epoch moment.
    ``model.jdsatepoch``
        Julian date of the epoch (computed from ``epochyr`` and ``epochdays``).
    ``model.ndot``
        First time derivative of the mean motion (ignored by SGP4).
    ``model.nddot``
        Second time derivative of the mean motion (ignored by SGP4).
    ``model.bstar``
        Ballistic drag coefficient B* in inverse earth radii.
    ``model.inclo``
        Inclination in radians.
    ``model.nodeo``
        Right ascension of ascending node in radians.
    ``model.ecco``
        Eccentricity.
    ``model.argpo``
        Argument of perigee in radians.
    ``model.mo``
        Mean anomaly in radians.
    ``model.no``
        Mean motion in radians per minute.

    """
    center = 399
    center_name = '399 EARTH'

    def __init__(self, line1, line2, name=None, ts=None):
        ts = ts or _ts

        self.name = None if name is None else name.strip()
        sat = Satrec.twoline2rv(line1, line2)
        self.model = sat

        # TODO: just use the Julian dates instead
        two_digit_year = sat.epochyr
        if two_digit_year < 57:
            year = two_digit_year + 2000;
        else:
            year = two_digit_year + 1900;

        self.epoch = ts.utc(year, 1, sat.epochdays)

        self.target = -100000 - self.model.satnum
        self.target_name = 'Satellite{0} {1}'.format(
            self.model.satnum,
            ' ' + repr(self.name) if self.name else '',
        )

    def __str__(self):
        sat = self.model
        return 'EarthSatellite{0} number={1!r} epoch={2}'.format(
            ' ' + repr(self.name) if self.name else '',
            sat.satnum,
            self.epoch.utc_iso(),
        )

    def __repr__(self):
        return '<{0}>'.format(self)

    def _position_and_velocity_TEME_km(self, t):
        """Return the raw true equator mean equinox (TEME) vectors from SGP4.

        Returns a tuple of NumPy arrays ``([x y z], [xdot ydot zdot])``
        expressed in kilometers and kilometers per second.  Note that we
        assume the TLE epoch to be a UTC date, per AIAA 2006-6753.

        """
        sat = self.model
        jd = t._utc_float()
        if getattr(jd, 'shape', None):
            e, r, v = sat.sgp4_array(jd, zeros_like(jd))
            messages = [SGP4_ERRORS[error] if error else None for error in e]
            return r.T, v.T, messages
        else:
            error, position, velocity = sat.sgp4(jd, 0.0)
            message = SGP4_ERRORS[error] if error else None
            return array(position), array(velocity), message

    def ITRF_position_velocity_error(self, t):
        """Return the ITRF position, velocity, and error at time `t`.

        The position is an x,y,z vector measured in au, the velocity is
        an x,y,z vector measured in au/day, and the error is a vector of
        possible error messages for the time or vector of times `t`.

        """
        rTEME, vTEME, error = self._position_and_velocity_TEME_km(t)
        rTEME /= AU_KM
        vTEME /= AU_KM
        vTEME *= DAY_S
        rITRF, vITRF = TEME_to_ITRF(t.ut1, rTEME, vTEME)
        return rITRF, vITRF, error

    def _at(self, t):
        """Compute this satellite's GCRS position and velocity at time `t`."""
        rITRF, vITRF, error = self.ITRF_position_velocity_error(t)
        rGCRS, vGCRS = ITRF_to_GCRS2(t, rITRF, vITRF)
        return rGCRS, vGCRS, rGCRS, error

    def find_events(self, topos, t0, t1, minimum_altitude_degrees=0.0):
        """Return the times at which the satellite rises, culminates, and sets.

        Searches between ``t0`` and ``t1``, which should each be a
        Skyfield :class:`~skyfield.timelib.Time` object, for passes of
        this satellite above the location ``topos`` that reach at least
        ``minimum_altitude_degrees`` above the horizon.

        Returns a tuple ``(t, events)`` whose first element is a
        :class:`~skyfield.timelib.Time` array and whose second element
        is an array of events:

        * 0 — Satellite rose above ``minimum_altitude_degrees``.
        * 1 — Satellite culminated and started to descend again.
        * 2 — Satellite fell below ``minimum_altitude_degrees``.

        Note that multiple culminations in a row are possible when,
        without setting, the satellite reaches a second peak altitude
        after descending partway down the sky from the first one.

        """
        # First, we find the moments of maximum altitude over the time
        # period.  Some of these maxima will be negative, meaning the
        # satellite failed to crest the horizon.

        ts = t0.ts
        at = (self - topos).at
        half_second = 0.5 / DAY_S
        orbits_per_minute = self.model.no / tau
        orbits_per_day = 24 * 60 * orbits_per_minute
        rough_period = 1 / orbits_per_day

        # Long-period satellites might rise each day not because of
        # their own motion, but because the Earth rotates under them, so
        # check position at least each quarter-day.  We might need to
        # tighten this even further if experience someday shows it
        # missing a pass of a particular satellite.
        if rough_period > 0.25:
            rough_period = 0.25

        def cheat(t):
            """Avoid computing expensive values that cancel out anyway."""
            t.gast = t.tt * 0.0
            t.M = t.MT = _identity

        def altitude_at(t):
            cheat(t)
            return at(t).altaz()[0].degrees

        altitude_at.rough_period = rough_period
        tmax, altitude = _find_maxima(t0, t1, altitude_at, half_second, 12)
        if not tmax:
            return tmax, ones_like(tmax)

        # Next, filter out the maxima that are not high enough.

        keepers = altitude >= minimum_altitude_degrees
        jdmax = tmax.tt[keepers]
        ones = ones_like(jdmax, 'uint8')

        # Finally, find the rising and setting that bracket each maximum
        # altitude.  We guess that the satellite will be back below the
        # horizon in between each pair of adjancent maxima.

        def below_horizon_at(t):
            cheat(t)
            return at(t).altaz()[0].degrees < minimum_altitude_degrees

        # The `jdo` array are the times of maxima, with their averages
        # in between them.  The start and end times are thrown in too,
        # in case a rising or setting is lingering out between a maxima
        # and the ends of our range.  Could this perhaps still miss a
        # stubborn rising or setting near the ends?
        doublets = repeat(concatenate(((t0.tt,), jdmax, (t1.tt,))), 2)
        jdo = (doublets[:-1] + doublets[1:]) / 2.0

        trs, rs = _find_discrete(t0.ts, jdo, below_horizon_at, half_second, 8)

        jd = concatenate((jdmax, trs.tt))
        v = concatenate((ones, rs * 2))

        i = jd.argsort()
        return ts.tt_jd(jd[i]), v[i]

_second = 1.0 / (24.0 * 60.0 * 60.0)

def theta_GMST1982(jd_ut1):
    """Return the angle of Greenwich Mean Standard Time 1982 given the JD.

    This angle defines the difference between the idiosyncratic True
    Equator Mean Equinox (TEME) frame of reference used by SGP4 and the
    more standard Pseudo Earth Fixed (PEF) frame of reference.

    From AIAA 2006-6753 Appendix C.

    """
    t = (jd_ut1 - T0) / 36525.0
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t
    dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t
    theta = (jd_ut1 % 1.0 + g * _second % 1.0) * tau
    theta_dot = (1.0 + dg * _second / 36525.0) * tau
    return theta, theta_dot

def TEME_to_ITRF(jd_ut1, rTEME, vTEME, xp=0.0, yp=0.0):
    """Convert TEME position and velocity into standard ITRS coordinates.

    This converts a position and velocity vector in the idiosyncratic
    True Equator Mean Equinox (TEME) frame of reference used by the SGP4
    theory into vectors into the more standard ITRS frame of reference.
    The velocity should be provided in units per day, not per second.

    From AIAA 2006-6753 Appendix C.

    """
    theta, theta_dot = theta_GMST1982(jd_ut1)
    zero = theta_dot * 0.0
    angular_velocity = array([zero, zero, -theta_dot])
    R = rot_z(-theta)

    if len(rTEME.shape) == 1:
        rPEF = (R).dot(rTEME)
        vPEF = (R).dot(vTEME) + cross(angular_velocity, rPEF)
    else:
        rPEF = einsum('ij...,j...->i...', R, rTEME)
        vPEF = einsum('ij...,j...->i...', R, vTEME) + cross(
            angular_velocity, rPEF, 0, 0).T

    if xp == 0.0 and yp == 0.0:
        rITRF = rPEF
        vITRF = vPEF
    else:
        W = (rot_x(yp)).dot(rot_y(xp))
        rITRF = (W).dot(rPEF)
        vITRF = (W).dot(vPEF)
    return rITRF, vITRF
