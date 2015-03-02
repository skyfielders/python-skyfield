"""An interface between JPL ephemerides and Skyfield."""

import jplephem
from jplephem.spk import SPK
from jplephem.names import target_names
from numpy import max, min

from .chaining import Body, Segment
from .constants import AU_KM, C_AUDAY, DAY_S
from .functions import length_of
from .positionlib import Astrometric, Barycentric, Topos
from .timelib import takes_julian_date


class Kernel(dict):
    def __init__(self, file):
        if isinstance(file, str):
            file = open(file, 'rb')
        self.spk = SPK(file)

        segments = [Segment(s.center, s.target, _build_compute(s))
                    for s in self.spk.segments]
        codes = set(s.center for s in segments).union(
                    s.target for s in segments)

        for code in codes:
            body = Body(code, segments)
            self[code] = body
            raw_name = target_names.get(code, None)
            if raw_name is None:
                continue
            name = raw_name.lower().replace(' ', '_')
            setattr(self, name, body)


def _build_compute(segment):
    """Build a Skyfield `compute` callback for the SPK `segment`."""

    if segment.data_type == 2:
        def compute(jd):
            position, velocity = segment.compute_and_differentiate(jd.tdb)
            return position / AU_KM, velocity / AU_KM

    elif segment.data_type == 3:
        def compute(jd):
            six = segment.compute(jd.tdb)
            return six[:3] / AU_KM, six[3:] * DAY_S / AU_KM

    else:
        raise ValueError('SPK data type {} not yet supported segment'
                         .format(segment.data_type))
    return compute


# The older ephemerides that the code below tackles use a different
# value for the AU, so, for now (until we fix our tests?):

class Planet(object):
    def __init__(self, ephemeris, jplephemeris, jplname):
        self.ephemeris = ephemeris
        self.jplephemeris = jplephemeris
        self.jplname = jplname

    def __repr__(self):
        return '<Planet %s>' % (self.jplname,)

    @takes_julian_date
    def __call__(self, jd):
        """Return the x,y,z position of this planet at the given time."""
        position, velocity = self._position_and_velocity(jd.tdb)
        i = Barycentric(position, velocity, jd)
        i.ephemeris = self.ephemeris
        return i

    def _position(self, jd_tdb):
        e = self.jplephemeris
        c = e.position
        if self.jplname == 'earth':
            p = c('earthmoon', jd_tdb) - c('moon', jd_tdb) * e.earth_share
        elif self.jplname == 'moon':
            p = c('earthmoon', jd_tdb) + c('moon', jd_tdb) * e.moon_share
        else:
            p = c(self.jplname, jd_tdb)
        p /= AU_KM
        if getattr(jd_tdb, 'shape', ()) == ():
            # Skyfield, unlike jplephem, is willing to accept and return
            # plain scalars instead of only trafficking in NumPy arrays.
            p = p[:,0]
        return p

    def _position_and_velocity(self, jd_tdb):
        e = self.jplephemeris
        c = e.compute
        if self.jplname == 'earth':
            pv = c('earthmoon', jd_tdb) - c('moon', jd_tdb) * e.earth_share
        elif self.jplname == 'moon':
            pv = c('earthmoon', jd_tdb) + c('moon', jd_tdb) * e.moon_share
        else:
            pv = c(self.jplname, jd_tdb)
        pv /= AU_KM
        if getattr(jd_tdb, 'shape', ()) == ():
            # Skyfield, unlike jplephem, is willing to accept and return
            # plain scalars instead of only trafficking in NumPy arrays.
            pv = pv[:,0]
        return pv[:3], pv[3:]

    def _observe_from_bcrs(self, observer):
        # TODO: should also accept another ICRS?

        jd_tdb = observer.jd.tdb
        lighttime0 = 0.0
        position, velocity = self._position_and_velocity(jd_tdb)
        vector = position - observer.position.au
        euclidian_distance = distance = length_of(vector)

        for i in range(10):
            lighttime = distance / C_AUDAY
            delta = lighttime - lighttime0
            if -1e-12 < min(delta) and max(delta) < 1e-12:
                break
            lighttime0 = lighttime
            position, velocity = self._position_and_velocity(jd_tdb - lighttime)
            vector = position - observer.position.au
            distance = length_of(vector)
        else:
            raise ValueError('observe_from() light-travel time'
                             ' failed to converge')

        g = Astrometric(vector, velocity - observer.velocity.au_per_d,
                        observer.jd)
        g.observer = observer
        g.distance = euclidian_distance
        g.lighttime = lighttime
        return g

class Earth(Planet):

    def topos(self, latitude=None, longitude=None, latitude_degrees=None,
              longitude_degrees=None, elevation_m=0.0):
        """Return a ``Topos`` object for a specific location on Earth."""
        t = Topos(latitude, longitude, latitude_degrees,
                  longitude_degrees, elevation_m)
        t.ephemeris = self.ephemeris
        return t

    def satellite(self, text):
        from .sgp4lib import EarthSatellite
        lines = text.splitlines()
        return EarthSatellite(lines, self)

class Ephemeris(object):

    def __init__(self, module):

        self.jplephemeris = jplephem.Ephemeris(module)

        self.sun = Planet(self, self.jplephemeris, 'sun')
        self.mercury = Planet(self, self.jplephemeris, 'mercury')
        self.venus = Planet(self, self.jplephemeris, 'venus')
        self.earth = Earth(self, self.jplephemeris, 'earth')
        self.moon = Planet(self, self.jplephemeris, 'moon')
        self.mars = Planet(self, self.jplephemeris, 'mars')
        self.jupiter = Planet(self, self.jplephemeris, 'jupiter')
        self.saturn = Planet(self, self.jplephemeris, 'saturn')
        self.uranus = Planet(self, self.jplephemeris, 'uranus')
        self.neptune = Planet(self, self.jplephemeris, 'neptune')
        self.pluto = Planet(self, self.jplephemeris, 'pluto')

    def _position(self, name, jd):
        return getattr(self, name)._position(jd)

    def _position_and_velocity(self, name, jd):
        return getattr(self, name)._position_and_velocity(jd)
