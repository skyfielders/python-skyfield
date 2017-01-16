"""Vector functions and their composition."""

from numpy import max, min
from .constants import C_AUDAY
from .errors import raise_error_for_deprecated_time_arguments
from .functions import length_of
from .positionlib import build_position
from .timelib import Time

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
        if not isinstance(t, Time):
            raise ValueError('please provide the at() method with a Time'
                             ' instance as its argument, instead of the'
                             ' value {0!r}'.format(t))
        observer_data = ObserverData()
        p, v, observer_data.gcrs_position = self._at(t)
        self._snag_observer_data(observer_data, t)
        center = self.center
        if center == 0:
            observer_data.bcrs_position = p
            observer_data.bcrs_velocity = v
        position = build_position(p, v, t, center, self.target, observer_data)
        return position

    def _snag_observer_data(self, data, t):
        pass

    def _observe_from_bcrs(self, observer):
        assert self.center == 0
        return _correct_for_light_travel_time(observer, self)

    def geometry_of(self, other):
        # TODO: deprecate this
        if isinstance(other, str):
            # TODO: is this always the right ephemeris to use?
            other = self.positives[0].ephemeris[other]
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
        p, v = 0.0, 0.0
        for segment in self.positives:
            p2, v2, gcrs_position = segment._at(t)
            p += p2
            v += v2
        for segment in self.negatives:
            p2, v2, gcrs_position = segment._at(t)
            p -= p2
            v -= v2
        return p, v, gcrs_position

    def _snag_observer_data(self, observer_data, t):
        # TODO: does it make sense to go through both all the positives
        # and all the negatives?
        for segment in self.positives:
            segment._snag_observer_data(observer_data, t)
        for segment in self.negatives:
            segment._snag_observer_data(observer_data, t)


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
    tposition, tvelocity, gcrs_position = target._at(t)
    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    t_tdb = t.tdb
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if -1e-12 < min(delta) and max(delta) < 1e-12:
            break
        t2 = ts.tdb(jd=t_tdb - light_time)
        tposition, tvelocity, gcrs_position = target._at(t2)
        distance = length_of(tposition - cposition)
        light_time0 = light_time
    else:
        raise ValueError('light-travel time failed to converge')
    return tposition - cposition, tvelocity - cvelocity, light_time


class ObserverData(object):
    """Extra information about an observer."""
    # TODO: expand the documentation for this class

    __slots__ = ('altaz_rotation', 'elevation_m', 'ephemeris',
                 'gcrs_position', 'bcrs_position', 'bcrs_velocity')

    def __init__(self):
        self.altaz_rotation = None  #go ahead and precompute in case needed N
        self.elevation_m = None  #just keep segment then just keep observer?
        self.ephemeris = None  #keep this on observer instead?
        self.gcrs_position = None
        self.bcrs_position = None
        self.bcrs_velocity = None
