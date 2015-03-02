"""Compute position by chaining together answers from different ephemerides."""

from collections import namedtuple
from numpy import max, min

from .constants import C_AUDAY
from .functions import length_of
from .positionlib import Barycentric, ICRS
from .timelib import JulianDate, takes_julian_date

Segment = namedtuple('Segment', 'center target compute')


class Body(object):
    def __init__(self, code, segments):
        self.code = code
        self.segments = segments

    def geometry_of(self, body):
        center_chain, target_chain = _connect(self, body)
        return Geometry(self.code, body.code, center_chain, target_chain)

    def observe(self, body):
        every = self.segments + body.segments
        segment_dict = {segment.target: segment for segment in every}
        center_chain = list(_center(self.code, segment_dict))[::-1]
        target_chain = list(_center(body.code, segment_dict))[::-1]
        if not center_chain[0].center == target_chain[0].center == 0:
            raise ValueError('cannot observe() unless both bodies can be'
                             ' referenced to the solar system barycenter')
        return Solution(self.code, body.code, center_chain, target_chain)


def _connect(body1, body2):
    """Return ``(sign, segment)`` tuple list leading from body1 to body2."""
    every = body1.segments + body2.segments
    segment_dict = {segment.target: segment for segment in every}
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

    @takes_julian_date
    def at(self, jd):
        """Return the geometric Cartesian position and velocity."""
        pos, vel = _tally(self.center_chain, self.target_chain, jd)
        cls = Barycentric if self.center == 0 else ICRS
        return cls(pos, vel, jd)


class Solution(object):
    def __init__(self, center, target, center_chain, target_chain):
        self.center = center
        self.target = target
        self.center_chain = center_chain
        self.target_chain = target_chain

    @takes_julian_date
    def at(self, jd):
        """Return a light-time corrected astrometric position and velocity."""
        cposition, cvelocity = _tally([], self.center_chain, jd)
        tposition, tvelocity = _tally([], self.target_chain, jd)
        distance = length_of(tposition - cposition)
        lighttime0 = 0.0
        jd_tdb = jd.tdb
        for i in range(10):
            lighttime = distance / C_AUDAY
            delta = lighttime - lighttime0
            if -1e-12 < min(delta) and max(delta) < 1e-12:
                break
            jd2 = JulianDate(tdb=jd_tdb - lighttime)
            tposition, tvelocity = _tally([], self.target_chain, jd2)
            distance = length_of(tposition - cposition)
            lighttime0 = lighttime
        else:
            raise ValueError('observe_from() light-travel time'
                             ' failed to converge')
        cls = Barycentric if self.center == 0 else ICRS
        return cls(tposition - cposition, tvelocity - cvelocity, jd)


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
