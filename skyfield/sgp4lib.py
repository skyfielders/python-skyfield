"""An object representing an Earth-orbiting satellite."""

from numpy import array, cross
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

from .constants import C_AUDAY, DAY_S, KM_AU, T0, tau
from .functions import length_of, rot_x, rot_y, rot_z
from .positionlib import Apparent, ICRS, ITRF_to_GCRS
from .timescales import takes_julian_date

_minutes_per_day = 1440.

class EarthSatellite(object):
    """An Earth satellite loaded from a TLE file and propagated with SGP4."""

    def __init__(self, lines, earth):
        self.earth = earth
        self.sgp4_satellite = twoline2rv(*lines[-2:], whichconst=wgs72)

    def _position_and_velocity_TEME_km(self, jd):
        """Return the raw true equator mean equinox (TEME) vectors from SGP4.

        Returns a tuple of NumPy arrays ``([x y z], [xdot ydot zdot])``
        expressed in kilometers and kilometers per second.

        """
        minutes_past_epoch = (jd.ut1 - self.sgp4_satellite.jdsatepoch) * 1440.
        position, velocity = sgp4(self.sgp4_satellite, minutes_past_epoch)
        return (array(position), array(velocity))

    def _compute_GCRS(self, jd):
        """Compute where satellite is in space on a given date."""

        rTEME, vTEME = self._position_and_velocity_TEME_km(jd)
        rTEME *= KM_AU
        vTEME *= KM_AU * DAY_S

        rITRF, vITRF = TEME_to_ITRF(jd.ut1, rTEME, vTEME)
        rGCRS = ITRF_to_GCRS(jd, rITRF)
        vGCRS = array((0.0, 0.0, 0.0))  # todo: someday also compute vGCRS?

        return rGCRS, vGCRS

    def observe_from(self, observer):
        # TODO: what if someone on Mars tries to look at the ISS?

        rGCRS, vGCRS = self._compute_GCRS(observer.jd)
        g = Apparent(rGCRS - observer.rGCRS,
                     vGCRS - observer.vGCRS,
                     observer.jd)
        g.observer = observer
        # g.distance = euclidian_distance
        return g


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
    angular_velocity = array([0, 0, -theta_dot])
    R = rot_z(-theta)
    rPEF = (R).dot(rTEME)
    vPEF = (R).dot(vTEME) + cross(angular_velocity, rPEF)
    if xp == 0.0 and yp == 0.0:
        rITRF = rPEF
        vITRF = vPEF
    else:
        W = (rot_x(-yp)).dot(rot_y(-xp))
        rITRF = (W).dot(rPEF)
        vITRF = (W).dot(vPEF)
    return rITRF, vITRF
