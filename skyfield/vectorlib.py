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
                # TODO: put an explanation here
                raise ValueError()

        selfp = getattr(self, 'positives', None) or (self,)
        selfn = getattr(self, 'negatives', ())

        otherp = getattr(other, 'positives', None) or (other,)
        othern = getattr(other, 'negatives', ())

        return VectorSum(self.center, other.target,
                         self.center_name, other.target_name,
                         selfp + otherp, selfn + othern)

    def __sub__(self, other):
        if self.center != other.center:
            # TODO: put an explanation here
            raise ValueError()

        selfp = getattr(self, 'positives', None) or (self,)
        selfn = getattr(self, 'negatives', ())

        otherp = getattr(other, 'positives', None) or (other,)
        othern = getattr(other, 'negatives', ())

        return VectorSum(self.target, other.target,
                         self.target_name, other.target_name,
                         selfp + othern, selfn + otherp)

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        p, v = self._at(t)
        data = ObserverData()
        self._snag_observer_data(data, t)
        position = build_position(p, v, t, self.center, self.target, data)
        return position

    def _snag_observer_data(self, data, t):
        pass

    def _observe_from_bcrs(self, observer):
        assert self.center == 0
        return _correct_for_light_travel_time(observer, self)

    def geometry_of(self, other):
        # TODO: deprecate this
        if isinstance(other, str):
            other = self.first.ephemeris[other]
        return other - self

    def topos(self, latitude=None, longitude=None, latitude_degrees=None,
              longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):
        # TODO: deprecate this
        from .toposlib import Topos
        t = Topos(latitude, longitude, latitude_degrees,
                  longitude_degrees, elevation_m, x, y)
        return self + t

    def satellite(self, text):
        # TODO: deprecate this
        from .sgp4lib import EarthSatellite
        sat = EarthSatellite(text.splitlines())
        return self + sat


class VectorSum(VectorFunction):
    def __init__(self, center, target, center_name, target_name,
                 positives, negatives):
        self.center = center
        self.target = target
        self.center_name = center_name
        self.target_name = target_name
        self.positives = positives
        self.negatives = negatives
        self.first = positives[0]
        self.rest = positives[1:]
        self.ephemeris = getattr(positives[0], 'ephemeris', None)

    def __str__(self):
        positives = self.positives
        negatives = self.negatives
        lines = [' + ' + str(segment) for segment in positives]
        lines.extend(' - ' + str(segment) for segment in negatives)
        return 'Sum of {0} vectors:\n{1}'.format(
            len(positives) + len(negatives),
            '\n'.join(lines),
        )

    def __repr__(self):
        return '<{0} of {1} vectors {2} -> {3}>'.format(
            type(self).__name__,
            len(self.positives) + len(self.negatives),
            self.center_name,
            self.target_name,
        )

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
        # TODO: does it make sense to go through both all the positives
        # and all the negatives?  This is trying to cover both the case
        # my_topos.at(t) where the topos is the final positive segment,
        # and also the case (satellite - my_topos).at(t) where the topos
        # is the negative segment.
        # TODO: fix this crazy gcrs business, which is computing again a
        # quantity that we already computed while doing the sum or diff.
        for segment in self.positives:
            segment._snag_observer_data(data, t)
            if segment.center == 399:
                p, v = segment._at(t)
                data.gcrs_position = p
        for segment in self.negatives:
            segment._snag_observer_data(data, t)
            if segment.center == 399:
                p, v = segment._at(t)
                data.gcrs_position = p


def _correct_for_light_travel_time(observer, target):
    """Return a light-time corrected astrometric position and velocity.

    Given an `observer` that is a `Barycentric` position somewhere in
    the solar system, compute where in the sky they will see the body
    `target`, by computing the light-time between them and figuring out
    where `target` was back when the light was leaving it that is now
    reaching the eyes or instruments of the `observer`.

    """
    t = observer.t
    ts = t.ts
    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d
    # TODO: .at() is a much more expensive operation than we really need
    # here; pivot to using ._at() here and in the second call down in the loop.
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
        raise ValueError('light-travel time failed to converge')
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
