"""Vector functions and their composition."""

from numpy import max, min
from .constants import C_AUDAY
from .errors import raise_error_for_deprecated_time_arguments
from .functions import length_of
from .positionlib import Astrometric, build_position

class VectorFunction(object):
    def __add__(self, other):
        if self.target != other.center:
            if other.target == self.center:
                self, other = other, self
            else:
                raise ValueError()

        selfp = getattr(self, 'positives', None) or (self,)
        selfn = getattr(self, 'negatives', ())

        otherp = getattr(other, 'positives', None) or (other,)
        othern = getattr(other, 'negatives', ())

        return Sum(self.center, other.target, selfp + otherp, selfn + othern)

    def __sub__(self, other):
        if self.center != other.center:
            raise ValueError()

        selfp = getattr(self, 'positives', None) or (self,)
        selfn = getattr(self, 'negatives', ())

        otherp = getattr(other, 'positives', None) or (other,)
        othern = getattr(other, 'negatives', ())

        return Sum(self.target, other.target, selfp + othern, selfn + otherp)

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        p, v = self._at(t)
        position = build_position(p, v, t, self.center, self.target)
        position.observer_data = data = ObserverData()
        self._snag_observer_data(data, t)
        return position

    def _snag_observer_data(self, data, t):
        pass

    def _observe_from_bcrs(self, observer):
        assert self.center == 0
        return observe(observer, self)

    def geometry_of(self, other):
        # TODO: deprecate this
        if isinstance(other, str):
            other = self.first.ephemeris.get(other)
        return other - self

    def topos(self, latitude=None, longitude=None, latitude_degrees=None,
              longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):
        # TODO: deprecate this
        assert self.target == 399
        from .toposlib import Topos
        t = Topos(latitude, longitude, latitude_degrees,
                  longitude_degrees, elevation_m, x, y)
        return self + t
        # t.ephemeris = self.ephemeris
        # t.segments += self.segments
        # return t


class Sum(VectorFunction):
    def __init__(self, center, target, positives, negatives):
        self.center = center
        self.target = target
        self.positives = positives
        self.negatives = negatives
        self.first = positives[0]
        self.rest = positives[1:]
        self.ephemeris = getattr(positives[0], 'ephemeris')

    def __str__(self):
        positives = self.positives
        negatives = self.negatives
        lines = [' + ' + str(segment) for segment in positives]
        lines.extend(' - ' + str(segment) for segment in negatives)
        return 'Sum of {} vectors:\n{}'.format(
            len(positives) + len(negatives),
            '\n'.join(lines),
        )

    def __repr__(self):
        segments = self.positives + self.negatives
        return '<Sum of {}>'.format(' '.join(repr(s) for s in segments))

    def _at(self, t):
        p, v = self.first._at(t)
        for segment in self.rest:
            p2, v2 = segment._at(t)
            p += p2
            v += v2
        for segment in self.negatives:
            p2, v2 = segment._at(t)
            p -= p2
            v -= v2
        return p, v

    def _snag_observer_data(self, data, t):
        for segment in self.positives:
            segment._snag_observer_data(data, t)
        # TODO: fix the fact we are having to compute this twice. Ugh.
        if segment.center == 399:
            p, v = segment._at(t)
            data.gcrs_position = p



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
    t = observer.t
    ts = t.ts
    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d
    t_bary = target.at(t)
    tposition = t_bary.position.au
    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    t_tdb = t.tdb
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if -1e-12 < min(delta) and max(delta) < 1e-12:
            break
        t2 = ts.tdb(jd=t_tdb - light_time)
        t_bary = target.at(t2)
        tposition = t_bary.position.au
        distance = length_of(tposition - cposition)
        light_time0 = light_time
    else:
        raise ValueError('observe() light-travel time failed to converge')
    tvelocity = t_bary.velocity.au_per_d
    pos = Astrometric(tposition - cposition, tvelocity - cvelocity, t)
    pos.light_time = light_time
    pos.observer = observer
    return pos


class ObserverData(object):
    """Extra information about an observer."""
    # TODO: expand the documentation for this class

    def __init__(self):
        self.altaz_rotation = None
        self.elevation = None
        self.ephemeris = None
        self.gcrs_position = None
