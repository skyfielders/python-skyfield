"""Vector functions and their composition."""

from jplephem.names import target_names as _jpl_code_name_dict
from numpy import max
from .constants import C_AUDAY
from .errors import DeprecationError
from .functions import length_of
from .positionlib import build_position
from .timelib import Time

class VectorFunction(object):
    """Given a time, computes a corresponding position."""

    ephemeris = None

    @property
    def center_name(self):
        return _jpl_name(self.center)

    @property
    def target_name(self):
        return _jpl_name(self.target)

    def __add__(self, other):
        if self.target != other.center:
            if other.target == self.center:
                self, other = other, self
            else:
                raise ValueError(
                    "you can only add two vectors"
                    " if the target where one of the vectors ends"
                    " is the center where the other vector starts"
                )

        selfp = getattr(self, 'positives', None) or (self,)
        selfn = getattr(self, 'negatives', ())

        otherp = getattr(other, 'positives', None) or (other,)
        othern = getattr(other, 'negatives', ())

        return VectorSum(self.center, other.target,
                         selfp + otherp, selfn + othern)

    def __neg__(self):
        return ReversedVector(self)

    def __sub__(self, other):
        if self.center != other.center:
            raise ValueError(
                "you can only subtract two vectors"
                " if they both start at the same center"
            )

        selfp = getattr(self, 'positives', None) or (self,)
        selfn = getattr(self, 'negatives', ())

        otherp = getattr(other, 'positives', None) or (other,)
        othern = getattr(other, 'negatives', ())

        return VectorSum(other.target, self.target,
                         selfp + othern, selfn + otherp)

    def at(self, t):
        """At time ``t``, compute the target's position relative to the center.

        If ``t`` is an array of times, then the returned position object
        will specify as many positions as there were times.  The kind of
        position returned depends on the value of the ``center``
        attribute:

        * Solar System Barycenter: :class:`~skyfield.positionlib.Barycentric`
        * Center of the Earth: :class:`~skyfield.positionlib.Geocentric`
        * Difference: :class:`~skyfield.positionlib.Geometric`
        * Anything else: :class:`~skyfield.positionlib.ICRF`

        """
        if not isinstance(t, Time):
            raise ValueError('please provide the at() method with a Time'
                             ' instance as its argument, instead of the'
                             ' value {0!r}'.format(t))
        observer_data = ObserverData()
        p, v, observer_data.gcrs_position, message = self._at(t)
        center = self.center
        position = build_position(p, v, t, center, self.target, observer_data)
        position._ephemeris = self.ephemeris
        position.message = message
        return position

    def _observe_from_bcrs(self, observer):
        if self.center != 0:
            raise ValueError('you can only observe() a body whose vector'
                             ' center is the Solar System Barycenter,'
                             ' but this vector has the center {0}'
                             .format(self.center_name))
        return _correct_for_light_travel_time(observer, self)

    def geometry_of(self, other):
        raise DeprecationError(
"""the geometry_of() method has, alas, been deprecated

This old method has been replaced by an improved interface.  If you just
need your software working again, install Skyfield 0.9.1 for a quick fix:

    pip install skyfield==0.9.1

Or, to update your old code, replace each operation that looks like:

    position = boston.geometry_of(satellite).at(t)

with the vector math that was previously hiding inside the old method:

    position = (satellite - boston).at(t)""")

    def topos(self, latitude=None, longitude=None, latitude_degrees=None,
              longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):
        raise DeprecationError(
"""the topos() method has, alas, been deprecated

This old method has been replaced by an improved interface.  If you just
need your software working again, install Skyfield 0.9.1 for a quick fix:

    pip install skyfield==0.9.1

Or, to update your old code, replace each operation that looks like:

    boston = earth.topos(...)

with the vector math that was previously hiding inside the old method:

    from skyfield.api import Topos
    boston = earth + Topos(...)""")

    def satellite(self, text):
        raise DeprecationError(
"""the satellite() method has, alas, been deprecated

This old method has been replaced by an improved interface.  If you just
need your software working again, install Skyfield 0.9.1 for a quick fix:

    pip install skyfield==0.9.1

Or, to update your old code, replace each operation that looks like:

    sat = earth.satellite(tle_text)

with the vector math (and the little bit of text manipulation) that was
previously hiding inside the old method:

    from skyfield.api import EarthSatellite
    line1, line2 = tle_text.splitlines()[-2:]
    sat = earth + EarthSatellite(line1, line2)""")

class ReversedVector(VectorFunction):
    def __init__(self, vector_function):
        self.center = vector_function.target
        self.target = vector_function.center
        self.vector_function = vector_function

    def _at(self, t):
        p, v, gcrs_position, message = self.vector_function._at(t)
        return -p, -v, gcrs_position, message

class VectorSum(VectorFunction):
    def __init__(self, center, target, positives, negatives):
        self.center = center
        self.target = target
        self.positives = positives
        self.negatives = negatives

        # For now, just grab the first ephemeris we can find.
        ephemerides = (segment.ephemeris for segments in (positives, negatives)
                       for segment in segments if segment.ephemeris)
        self.ephemeris = next(ephemerides, None)

    def __str__(self):
        positives = self.positives
        negatives = self.negatives
        lines = [' - ' + str(segment) for segment in reversed(negatives)]
        lines.extend(' + ' + str(segment) for segment in positives)
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
            p2, v2, gcrs_position, message = segment._at(t)
            p += p2
            v += v2
        for segment in self.negatives:
            p2, v2, gcrs_position, ignored_message = segment._at(t)
            p -= p2
            v -= v2
        return p, v, gcrs_position, message

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
    whole = t.whole
    tdb_fraction = t.tdb_fraction

    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d

    tposition, tvelocity, gcrs_position, message = target._at(t)

    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if abs(max(delta)) < 1e-12:
            break

        # We assume a light travel time of at most a couple of days.  A
        # longer light travel time would best be split into a whole and
        # fraction, for adding to the whole and fraction of TDB.
        t2 = ts.tdb_jd(whole, tdb_fraction - light_time)

        tposition, tvelocity, gcrs_position, message = target._at(t2)
        distance = length_of(tposition - cposition)
        light_time0 = light_time
    else:
        raise ValueError('light-travel time failed to converge')
    return tposition - cposition, tvelocity - cvelocity, t, light_time


class ObserverData(object):
    """Essential facts about an observer, that may be needed later."""

    __slots__ = ('elevation_m', 'gcrs_position')

    def __init__(self):
        self.elevation_m = None
        self.gcrs_position = None

def _jpl_name(code_or_string):
    if isinstance(code_or_string, int):
        name = _jpl_code_name_dict.get(code_or_string)
        if name is not None:
            return '{0} {1}'.format(code_or_string, name)
    return str(code_or_string)
