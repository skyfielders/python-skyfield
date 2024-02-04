"""Vector functions and their composition."""

from jplephem.names import target_names as _jpl_code_name_dict
from numpy import max
from .constants import C_AUDAY
from .descriptorlib import reify
from .errors import DeprecationError
from .functions import length_of
from .positionlib import build_position
from .timelib import Time

class VectorFunction(object):
    """Given a time, computes a corresponding position."""

    ephemeris = None

    @reify
    def vector_name(self):
        return type(self).__name__

    @reify
    def center_name(self):
        return _jpl_name(self.center)

    @reify
    def target_name(self):
        return _jpl_name(self.target)

    def __repr__(self):
        return '<{0} {1}>'.format(type(self).__name__, str(self))

    def __str__(self):
        if self.target is self:
            return self.target_name
        return self.arrow_str()

    def arrow_str(self):
        return '{0} {1} -> {2}'.format(
            self.vector_name, self.center_name, self.target_name,
        )

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

        self_vfs = getattr(self, 'vector_functions', None) or (self,)
        other_vfs = getattr(other, 'vector_functions', None) or (other,)

        return VectorSum(self.center, other.target, self_vfs + other_vfs)

    def __neg__(self):
        return ReversedVector(self)

    def __sub__(self, other):
        if self.center != other.center:
            raise ValueError(
                "you can only subtract two vectors"
                " if they both start at the same center"
            )

        self_vfs = getattr(self, 'vector_functions', None) or (self,)
        other_vfs = getattr(other, 'vector_functions', None) or (other,)
        other_vfs = tuple(reversed([-vf for vf in other_vfs]))

        return VectorSum(other.target, self.target, other_vfs + self_vfs)

    def at(self, t):
        """At time ``t``, compute the target's position relative to the center.

        If ``t`` is an array of times, then the returned position object
        will specify as many positions as there were times.  The kind of
        position returned depends on the value of the ``center``
        attribute:

        * Solar System Barycenter: :class:`~skyfield.positionlib.Barycentric`
        * Center of the Earth: :class:`~skyfield.positionlib.Geocentric`
        * Anything else: :class:`~skyfield.positionlib.ICRF`

        """
        if not isinstance(t, Time):
            raise ValueError('please provide the at() method with a Time'
                             ' instance as its argument, instead of the'
                             ' value {0!r}'.format(t))
        p, v, gcrs_position, message = self._at(t)
        center = self.center
        position = build_position(p, v, t, center, self.target)
        position._ephemeris = self.ephemeris
        position._observer_gcrs_au = gcrs_position
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

    @reify
    def vector_name(self):
        return 'Reversed ' + self.vector_function.vector_name

    @reify
    def center_name(self):
        return self.vector_function.target_name

    @reify
    def target_name(self):
        return self.vector_function.center_name

    def __neg__(self):
        return self.vector_function

    def _at(self, t):
        p, v, _, message = self.vector_function._at(t)
        return -p, -v, None, message

class VectorSum(VectorFunction):
    def __init__(self, center, target, vector_functions):
        self.center = center
        self.target = target
        self.vector_functions = vector_functions

        # For now, just grab the first ephemeris we can find.
        ephemerides = (segment.ephemeris for segment in vector_functions
                       if segment.ephemeris)
        self.ephemeris = next(ephemerides, None)

    def __str__(self):
        vector_functions = self.vector_functions
        lines = [' ' + segment.arrow_str() for segment in vector_functions]
        return 'Sum of {0} vectors:\n{1}'.format(
            len(vector_functions),
            '\n'.join(lines),
        )

    def __repr__(self):
        return '<Vector{0}>'.format(self)

    def _at(self, t):
        p, v = 0.0, 0.0
        gcrs_position = None
        vfs = self.vector_functions
        for vf in vfs:
            p2, v2, _, message = vf._at(t)
            if vf.center == 399:
                gcrs_position = -p
            p += p2
            v += v2
        if vfs[0].center == 0 and vf.center == 399:
            gcrs_position = p2
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
        if max(abs(delta), initial=0.0) < 1e-12:
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

def _jpl_name(target):
    if not isinstance(target, int):
        return type(target).__name__
    name = _jpl_code_name_dict.get(target)
    if name is None:
        return str(target)
    return '{0} {1}'.format(target, name)
