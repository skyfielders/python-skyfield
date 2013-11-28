"""An object representing an Earth-orbiting satellite."""

from numpy import array
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4

from .constants import C_AUDAY, DAY_S, KM_AU
from .functions import length_of
from .positionlib import ICRS, Astrometric
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

    @takes_julian_date
    def __call__(self, jd):
        """Compute where satellite is in space on a given date."""

        position_teme, velocity_teme = self._position_and_velocity_TEME_km(jd)
        # TODO: real conversion from TEME to GCRS
        position = array(position_teme) * KM_AU
        velocity = array(velocity_teme) * KM_AU * DAY_S
        e = self.earth(jd)
        p = ICRS(position + e.position, velocity + e.velocity, jd=jd)
        # a.ephemeris = self.ephemeris
        return p

    def observe_from(self, observer):
        # TODO: should we add the kind of light-time back-dating that
        # planets use, in case, say, someone on Mars tries to look at
        # the ISS?

        my = self(observer.jd)
        g = Astrometric(my.position - observer.position,
                        my.velocity - observer.velocity,
                        observer.jd)
        g.observer = observer
        # g.distance = euclidian_distance
        g.lighttime = length_of(g.position) / C_AUDAY
        return g
