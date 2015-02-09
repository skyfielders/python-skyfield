"""An interface between JPL ephemerides and Skyfield."""

import jplephem
from collections import defaultdict, deque
from jplephem.spk import SPK
from jplephem.names import target_names
from numpy import max, min

from .constants import AU_KM, C_AUDAY
from .functions import length_of
from .positionlib import Barycentric, Astrometric, Topos
from .timelib import takes_julian_date


class Kernel(object):
    def __init__(self, open_file):
        self.spk = SPK(open_file)

        edges = defaultdict(list)
        for segment in self.spk.segments:
            for code in segment.center, segment.target:
                edges[code].append(segment)

        self.edges = edges
        self.codes = set(edges.keys())
        self.leaves = set(code for code in self.codes if len(edges[code]) == 1)

        for code in self.codes:
            name = target_names[code].lower().replace(' ', '_')
            setattr(self, name, Body(self, code))


class Body(object):
    def __init__(self, kernel, code):
        self.kernel = kernel
        self.code = code
        self.targets = {}

    def observe(self, body):
        if body in self.targets:
            return self.targets

        if self.kernel is not body.kernel:
            raise ValueError('cross-kernel positions not yet implemented')

        here, there = self.code, body.code

        if here == there:
            raise ValueError('a body cannot observe itself')

        # For efficiency, we pretend to have already visited leaf nodes
        # because they, by definition, cannot move us toward the target.

        paths = dict.fromkeys(self.kernel.leaves)
        paths.pop(there, None)
        paths[here] = ()

        # Standard breadth-first search.

        places = deque()
        places.append(here)
        while places:
            here = places.popleft()
            for segment in self.kernel.edges[here]:
                code = _other(segment, here)
                if code not in paths:
                    paths[code] = paths[here] + (segment,)
                    if code == there:
                        return Solution(paths[code])
                    places.append(code)

        raise ValueError('{0} cannot observe {1}'.format(self.code, body.code))


def _other(segment, code):
    """Return the other code besides `code` that a segment names."""
    return segment.center if (segment.target == code) else segment.target


class Solution(object):
    def __init__(self, path):
        self.path = path


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
        vector = position - observer.position.AU
        euclidian_distance = distance = length_of(vector)

        for i in range(10):
            lighttime = distance / C_AUDAY
            delta = lighttime - lighttime0
            if -1e-12 < min(delta) and max(delta) < 1e-12:
                break
            lighttime0 = lighttime
            position, velocity = self._position_and_velocity(jd_tdb - lighttime)
            vector = position - observer.position.AU
            distance = length_of(vector)
        else:
            raise ValueError('observe_from() light-travel time'
                             ' failed to converge')

        g = Astrometric(vector, velocity - observer.velocity.AU_per_d,
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
