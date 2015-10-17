"""n interface between JPL ephemerides and Skyfield."""

from collections import namedtuple
from numpy import max, min

from jplephem.spk import SPK
from jplephem.names import target_names as _names

from .constants import AU_KM, C_AUDAY, DAY_S
from .functions import length_of
from .positionlib import Astrometric, Barycentric, ICRS
from .timelib import JulianDate, takes_julian_date
from .units import Distance, Velocity

Segment = namedtuple('Segment', 'center target compute')
_targets = dict((name, target) for (target, name) in _names.items())


class SpiceKernel(object):
    def __init__(self, filename):
        self.filename = filename
        self.spk = SPK.open(filename)
        self.segments = [Segment(s.center, s.target, _build_compute(s))
                         for s in self.spk.segments]
        self.codes = set(s.center for s in self.segments).union(
                         s.target for s in self.segments)

    def __str__(self):
        return str(self.spk)

    def __getitem__(self, name):
        code = self.decode(name)
        return Body(self, code)

    def decode(self, name):
        if isinstance(name, int):
            return name
        name = name.upper()
        code = _targets.get(name)
        if code is None:
            raise KeyError('unknown SPICE target name {0!r}'.format(name))
        if code not in self.codes:
            names = ', '.join(_names[c] for c in self.codes)
            raise KeyError('kernel {0} is missing {1!r} -'
                           ' the targets it supports are: {2}'
                           .format(self.filename, name, names))
        return code


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


class Body(object):
    def __init__(self, ephemeris, code):
        self.ephemeris = ephemeris
        self.segments = ephemeris.segments
        self.code = code

    @takes_julian_date
    def at(self, jd):
        """Compute the Solar System position of this body at a given time."""
        segments = self.segments
        segment_dict = dict((segment.target, segment) for segment in segments)
        chain = list(_center(self.code, segment_dict))[::-1]
        pos, vel = _tally((), chain, jd)
        barycentric = Barycentric(pos, vel, jd)
        barycentric.ephemeris = self.ephemeris
        return barycentric

    def geometry_of(self, body):
        if not isinstance(body, Body):
            code = self.ephemeris.decode(body)
            body = Body(self.ephemeris, code)
        center_chain, target_chain = _connect(self, body)
        return Geometry(self.code, body.code, center_chain, target_chain)

    def _observe_from_bcrs(self, observer):
        return observe(observer, self)

    def topos(self, latitude=None, longitude=None, latitude_degrees=None,
              longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):
        assert self.code == 399
        from .toposlib import Topos
        t = Topos(latitude, longitude, latitude_degrees,
                  longitude_degrees, elevation_m, x, y)
        t.ephemeris = self.ephemeris
        t.segments += self.segments
        return t

    def satellite(self, text):
        assert self.code == 399
        from .sgp4lib import EarthSatellite
        lines = text.splitlines()
        return EarthSatellite(lines, self)


def observe(observer, target):
    """Return a light-time corrected astrometric position and velocity.

    Given an `observer` that is a `Barycentric` position somewhere in
    the solar system, compute where in the sky they will see the body
    `target`, by computing the light-time between them and figuring out
    where `target` was back when the light was leaving it that is now
    reaching the eyes or instruments of the `observer`.

    """
    # cposition, cvelocity = _tally([], self.center_chain, jd)
    # tposition, tvelocity = _tally([], self.target_chain, jd)
    jd = observer.jd
    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d
    tjd = target.at(jd)
    tposition = tjd.position.au
    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    jd_tdb = jd.tdb
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if -1e-12 < min(delta) and max(delta) < 1e-12:
            break
        jd2 = JulianDate(tdb=jd_tdb - light_time)
        tjd = target.at(jd2)
        tposition = tjd.position.au
        distance = length_of(tposition - cposition)
        light_time0 = light_time
    else:
        raise ValueError('observe_from() light-travel time'
                         ' failed to converge')
    tvelocity = tjd.velocity.au_per_d
    pos = Astrometric(tposition - cposition, tvelocity - cvelocity, jd)
    pos.light_time = light_time
    class Observer(object):
        pass
    pos.observer = Observer()
    pos.observer.position = Distance(cposition)
    pos.observer.velocity = Velocity(cvelocity)
    pos.observer.geocentric = False  # TODO
    pos.observer.ephemeris = target.ephemeris
    if observer.altaz_rotation is not None:
        pos.observer.altaz_rotation = observer.altaz_rotation
        pos.observer.topos = observer.topos
    return pos


def _connect(body1, body2):
    """Return ``(sign, segment)`` tuple list leading from body1 to body2."""
    every = body1.segments + body2.segments
    segment_dict = dict((segment.target, segment) for segment in every)
    segments1 = list(_center(body1.code, segment_dict))[::-1]
    segments2 = list(_center(body2.code, segment_dict))[::-1]
    if segments1[0].center != segments2[0].center:
        raise ValueError('cannot trace these bodies back to a common center')
    i = sum(1 for s1, s2 in zip(segments1, segments2) if s1.target == s2.target)
    return segments1[i:], segments2[i:]


def _center(code, segment_dict):
    """Starting with `code`, follow segments from target to center."""
    while code in segment_dict:
        segment = segment_dict[code]
        yield segment
        code = segment.center


class Geometry(object):
    def __init__(self, center, target, center_chain, target_chain):
        self.center = center
        self.target = target
        self.center_chain = center_chain
        self.target_chain = target_chain

    def __str__(self):
        return 'Geometry\n{0}'.format('\n'.join(
            ' {0}'.format(c)
            for c in self.center_chain + self.target_chain))

    @takes_julian_date
    def at(self, jd):
        """Return the geometric Cartesian position and velocity."""
        pos, vel = _tally(self.center_chain, self.target_chain, jd)
        cls = Barycentric if self.center == 0 else ICRS
        return cls(pos, vel, jd)


def _tally(minus_chain, plus_chain, jd):
    position = velocity = 0.0
    for segment in minus_chain:
        p, v = segment.compute(jd)
        position -= p
        velocity -= v
    for segment in plus_chain:
        p, v = segment.compute(jd)
        position += p
        velocity += v
    return position, velocity
