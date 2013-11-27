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
    """A satellite in orbit around the Earth, initialized by a TLE entry."""

    def __init__(self, lines, earth):
        self.earth = earth
        self.sgp4_satellite = twoline2rv(*lines[1:3], whichconst=wgs72)

    @takes_julian_date
    def __call__(self, jd):
        """Compute where satellite is in space on a given date."""

        m = (jd.ut1 - self.sgp4_satellite.jdsatepoch) * _minutes_per_day
        position, velocity = sgp4(self.sgp4_satellite, m)
        position = array(position) * KM_AU
        velocity = array(velocity) * KM_AU * DAY_S
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
