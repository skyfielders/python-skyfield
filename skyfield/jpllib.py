"""An interface between JPL ephemerides and Skyfield."""

import jplephem
from collections import defaultdict, deque
from jplephem.spk import SPK
from jplephem.names import target_names
from numpy import max, min

from .constants import AU_KM, C_AUDAY, DAY_S
from .functions import length_of
from .positionlib import Astrometric, Barycentric, ICRS, Topos
from .timelib import takes_julian_date


class Kernel(dict):
    def __init__(self, file):
        if isinstance(file, str):
            file = open(file, 'rb')

        self.spk = SPK(file)

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

    def geometry_of(self, body):
        if body in self.targets:
            return self.targets

        if self.kernel is not body.kernel:
            raise ValueError('cross-kernel positions not yet implemented')

        path = _find_segments_connecting(self.kernel, self.code, body.code)
        chain = _build_chain(path, self.code)
        return Geometry(self.code, body.code, chain)

    def observe(self, body):
        if body in self.targets:
            return self.targets

        if self.kernel is not body.kernel:
            raise ValueError('cross-kernel positions not yet implemented')

        cpath = _find_segments_connecting(self.kernel, 0, self.code)
        tpath = _find_segments_connecting(self.kernel, 0, body.code)
        center_chain = _build_chain(cpath, 0)
        target_chain = _build_chain(tpath, 0)
        return Solution(self.code, body.code, center_chain, target_chain)


def _other(segment, code):
    """Return the other code besides `code` that a segment names."""
    return segment.center if (segment.target == code) else segment.target


class Geometry(object):
    def __init__(self, center, target, chain):
        self.center = center
        self.target = target
        self.chain = chain

    @takes_julian_date
    def at(self, jd):
        """Return the geometric cartesian position and velociy."""
        position, velocity = _tally_chain(self.chain, jd.tdb)
        cls = Barycentric if self.center == 0 else ICRS
        return cls(position, velocity, jd)


class Solution(object):
    def __init__(self, center, target, center_chain, target_chain):
        self.center = center
        self.target = target
        self.center_chain = center_chain
        self.target_chain = target_chain

    @takes_julian_date
    def at(self, jd):
        """Return a light-time corrected astrometric position and velocity."""
        jd_tdb = jd.tdb
        cposition, cvelocity = _tally_chain(self.center_chain, jd_tdb)
        tposition, tvelocity = _tally_chain(self.target_chain, jd_tdb)
        distance = length_of(tposition - cposition)
        lighttime0 = 0.0
        for i in range(10):
            lighttime = distance / C_AUDAY
            delta = lighttime - lighttime0
            if -1e-12 < min(delta) and max(delta) < 1e-12:
                break
            tposition, tvelocity = _tally_chain(self.target_chain,
                                                jd_tdb - lighttime)
            distance = length_of(tposition - cposition)
            lighttime0 = lighttime
        else:
            raise ValueError('observe_from() light-travel time'
                             ' failed to converge')
        cls = Barycentric if self.center == 0 else ICRS
        return cls(tposition - cposition, tvelocity - cvelocity, jd)


def _find_segments_connecting(kernel, center, target):
    """Return a path from `center` to `target` using the `kernel`."""

    here, there = center, target
    if here == there:
        raise ValueError('a body cannot observe itself')

    # For efficiency, we pretend to have already visited leaf nodes
    # because they, by definition, cannot move us toward the target.

    paths = dict.fromkeys(kernel.leaves)
    paths.pop(there, None)
    paths[here] = ()

    # Standard breadth-first search.

    places = deque()
    places.append(here)
    while places:
        here = places.popleft()
        for segment in kernel.edges[here]:
            code = _other(segment, here)
            if code not in paths:
                paths[code] = paths[here] + (segment,)
                if code == there:
                    return paths[there]
                places.append(code)

    raise ValueError('{0} cannot observe {1}'.format(center, target))


def _build_chain(path, center):
    """Return a chain of segments that should be added or subtracted."""
    chain = []
    for segment in path:
        if segment.center == center:
            chain.append((1.0, segment))
            center = segment.target
        else:
            chain.append((-1.0, segment))
            center = segment.center
    return chain


def _tally_chain(chain, jd_tdb):
    position = velocity = 0.0

    for sign, segment in chain:
        if segment.data_type == 2:
            p, v = segment.compute_and_differentiate(jd_tdb)
            position += sign * p
            velocity += sign * v
        elif segment.data_type == 3:
            six = sign * segment.compute(jd_tdb)
            position += six[:3]
            velocity += six[3:] * DAY_S
        else:
            raise ValueError('SPK data type {} not yet supported segment'
                             .format(segment.data_type))

    return position / AU_KM, velocity / AU_KM


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
