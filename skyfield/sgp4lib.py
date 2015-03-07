"""An interface between Skyfield and the Python ``sgp4`` library."""

from datetime import timedelta

from numpy import array, cross, einsum, pi, sqrt, zeros_like
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

from .constants import AU_KM, DAY_S, GE, T0, tau
from .earthlib import earth_radius_at_latitude
from .functions import rot_x, rot_y, rot_z
from .positionlib import Apparent, Geocentric, GCRS_to_Topos, ITRF_to_GCRS
from .timelib import JulianDate, takes_julian_date
from .units import Distance

# important ones:
# jdsatepoch
# bstar
# inclo - inclination
# nodeo - right ascension of ascending node
# ecco - eccentricity
# argpo - argument of perigee
# mo - mean anomaly
# no - mean motion

_minutes_per_day = 1440.

class EarthSatellite(object):
    """An Earth satellite loaded from a TLE file and propagated with SGP4."""

    def __init__(self, lines, earth):
        sat = twoline2rv(*lines[-2:], whichconst=wgs72)
        self._sgp4_satellite = sat
        self._earth = earth
        self.epoch = JulianDate(utc=(sat.epochyr, 1, sat.epochdays - 1.0))

    def __repr__(self):
        sat = self._sgp4_satellite
        return '<EarthSatellite number={1!r} epoch={0}>'.format(
            self.epoch.utc_iso(), sat.satnum)

    def _position_and_velocity_TEME_km(self, jd):
        """Return the raw true equator mean equinox (TEME) vectors from SGP4.

        Returns a tuple of NumPy arrays ``([x y z], [xdot ydot zdot])``
        expressed in kilometers and kilometers per second.  Note that we
        assume the TLE epoch to be a UTC date, per AIAA 2006-6753.

        """
        sat = self._sgp4_satellite
        epoch = sat.jdsatepoch
        minutes_past_epoch = (jd._utc_float() - epoch) * 1440.
        if getattr(minutes_past_epoch, 'shape', None):
            position = []
            velocity = []
            error = []
            for m in minutes_past_epoch:
                p, v = sgp4(sat, m)
                position.append(p)
                velocity.append(v)
                error.append(sat.error_message)
            return array(position).T, array(velocity).T, error
        else:
            position, velocity = sgp4(sat, minutes_past_epoch)
            return array(position), array(velocity), sat.error_message

    def _compute_GCRS(self, jd):
        """Compute where satellite is in space on a given date."""

        rTEME, vTEME, error = self._position_and_velocity_TEME_km(jd)
        rTEME /= AU_KM
        vTEME /= AU_KM
        vTEME *= DAY_S

        rITRF, vITRF = TEME_to_ITRF(jd.ut1, rTEME, vTEME)
        rGCRS = ITRF_to_GCRS(jd, rITRF)
        vGCRS = zeros_like(rGCRS)  # todo: someday also compute vGCRS?

        return rGCRS, vGCRS, error

    @takes_julian_date
    def elevation(self, jd):
        """Return elevation of this satellite over ground.

        Return value is Distance adjusted for approximate shape of the
        Earth.

        """
        latitude = self.over_location(jd).latitude
        earth_radius = earth_radius_at_latitude(latitude.radians)
        return Distance(m=self.gcrs(jd).distance().m - earth_radius)

    @takes_julian_date
    def gcrs(self, jd):
        """Return a GCRS position for this Earth satellite.

        Uses standard SGP4 theory to predict the satellite location.

        """
        position_au, velociy_au_per_d, error = self._compute_GCRS(jd)
        g = Geocentric(position_au, velociy_au_per_d, jd)
        g.sgp4_error = error
        return g

    @property
    def mean_motion(self):
        # revolutions per second
        return self._sgp4_satellite.no / 60

    @property
    def orbital_period(self):
        return timedelta(seconds=2 * pi / self.mean_motion)

    @takes_julian_date
    def orbital_speed(self, jd):
        """Return orbital speed of this satellite in meters per second.
        """
        elevation = self.gcrs(jd).distance().m
        semi_major_axis = self.semi_major_axis.m
        return sqrt(GE * ((2 / elevation) - (1 / semi_major_axis)))

    @takes_julian_date
    def over_location(self, jd):
        """Return a Topos instance for the point on the Earth over which
        this satellite will be at the given date.

        >>> from skyfield.api import earth, JulianDate
        >>> tle = ("ISS (ZARYA)\n"
        ... "1 25544U 98067A   15058.48161588  .00023857  00000-0  35618-3 0  9991\n"
        ... "2 25544  51.6478 271.9610 0008043  54.4849  18.6646 15.54887163930964")
        >>> iss = earth.satellite(tle)
        >>> topos = iss.over_location(JulianDate(utc=(2015, 2, 27, 22, 22, 0)))
        >>> topos.latitude
        <Angle 49deg 19' 48.3">
        >>> topos.longitude
        <Angle -155deg 37' 41.8">
        """
        return GCRS_to_Topos(self.gcrs(jd).position.km, jd)

    def _observe_from_bcrs(self, observer):
        # TODO: what if someone on Mars tries to look at the ISS?

        jd = observer.jd
        rGCRS, vGCRS, error = self._compute_GCRS(jd)
        rGCRS - observer.rGCRS
        vGCRS - observer.vGCRS
        g = Apparent(rGCRS - observer.rGCRS, vGCRS - observer.vGCRS, jd)
        g.sgp4_error = error
        g.observer = observer
        # g.distance = euclidian_distance
        return g

    @property
    def semi_major_axis(self):
        return Distance(m=(GE / (self.mean_motion ** 2)) ** (1/3))


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
