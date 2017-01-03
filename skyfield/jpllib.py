"""An interface between JPL ephemerides and Skyfield."""

import os
from collections import defaultdict
from numpy import max, min

from jplephem.spk import SPK
from jplephem.names import target_name_pairs, target_names as _names

from .constants import AU_KM, C_AUDAY, DAY_S
from .errors import DeprecationError, raise_error_for_deprecated_time_arguments
from .functions import length_of
from .positionlib import Astrometric, Barycentric, ICRF, build_position
from .timelib import calendar_date

_targets = dict((name, target) for (target, name) in target_name_pairs)

class SpiceKernel(object):
    """Ephemeris file in NASA .bsp format.

    A "Spacecraft and Planet Kernel" (SPK) file from NASA provides
    (x,y,z) coordinates for bodies in the Solar System like the Sun,
    planets, moons, and spacecraft.

    You can download a .bsp file yourself and use this class to open it,
    or use the Skyfield `load` function to automatically download a
    popular ephemeris.  Once loaded, you can print this object to the
    screen to see a report on the segments that it includes:

    >>> planets = load('de421.bsp')
    >>> print(planets)
    SPICE kernel file 'de421.bsp' has 15 segments
      JD 2414864.50 - JD 2471184.50  (1899-07-28 through 2053-10-08)
          0 -> 1    SOLAR SYSTEM BARYCENTER -> MERCURY BARYCENTER
          0 -> 2    SOLAR SYSTEM BARYCENTER -> VENUS BARYCENTER
          0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
          0 -> 4    SOLAR SYSTEM BARYCENTER -> MARS BARYCENTER
          0 -> 5    SOLAR SYSTEM BARYCENTER -> JUPITER BARYCENTER
          0 -> 6    SOLAR SYSTEM BARYCENTER -> SATURN BARYCENTER
          0 -> 7    SOLAR SYSTEM BARYCENTER -> URANUS BARYCENTER
          0 -> 8    SOLAR SYSTEM BARYCENTER -> NEPTUNE BARYCENTER
          0 -> 9    SOLAR SYSTEM BARYCENTER -> PLUTO BARYCENTER
          0 -> 10   SOLAR SYSTEM BARYCENTER -> SUN
          3 -> 301  EARTH BARYCENTER -> MOON
          3 -> 399  EARTH BARYCENTER -> EARTH
          1 -> 199  MERCURY BARYCENTER -> MERCURY
          2 -> 299  VENUS BARYCENTER -> VENUS
          4 -> 499  MARS BARYCENTER -> MARS

    To create a `Body` object for a target you are interested in, use
    square brackets and supply the target's name or integer code:

    >>> planets['earth']
    <Body 399 'EARTH' from kernel 'de421.bsp'>
    >>> planets[499]
    <Body 499 'MARS' from kernel 'de421.bsp'>

    """
    def __init__(self, path):
        self.path = path
        self.filename = os.path.basename(path)
        self.spk = SPK.open(path)
        self.segments = [Segment(self.filename, s) for s in self.spk.segments]
        self.codes = set(s.center for s in self.segments).union(
                         s.target for s in self.segments)

    def __repr__(self):
        return '<{} {!r}>'.format(type(self).__name__, self.path)

    def __str__(self):
        segments = self.spk.segments
        lines = ['SPICE kernel file {0!r} has {1} segments'
                 .format(self.filename, len(segments))]
        format_date = '{0}-{1:02}-{2:02}'.format
        start = end = None
        for s in segments:
            if start != s.start_jd or end != s.end_jd:
                start, end = s.start_jd, s.end_jd
                starts = format_date(*calendar_date(int(start)))
                ends = format_date(*calendar_date(int(end)))
                lines.append('  JD {0:.2f} - JD {1:.2f}  ({2} through {3})'
                             .format(start, end, starts, ends))
            lines.append(_format_segment(s))
        return '\n'.join(lines)

    def __getitem__(self, name):
        """Return a `Body` given a target name or integer."""
        code = self.decode(name)
        return Body(self, code)

    def comments(self):
        """Return the comments string of this kernel.

        The resulting string often contains embedded newlines, and is
        formatted for a human reader.

        >>> print(planets.comments())
        ; de421.bsp LOG FILE
        ;
        ; Created 2008-02-12/11:33:34.00.
        ...
        LEAPSECONDS_FILE    = naif0007.tls
        SPK_FILE            = de421.bsp
        ...

        """
        return self.spk.comments()

    def names(self):
        """Return all target names that are valid with this kernel.

        >>> pprint(planets.names())
        {0: ['SOLAR_SYSTEM_BARYCENTER', 'SSB', 'SOLAR SYSTEM BARYCENTER'],
         1: ['MERCURY_BARYCENTER', 'MERCURY BARYCENTER'],
         2: ['VENUS_BARYCENTER', 'VENUS BARYCENTER'],
         3: ['EARTH_BARYCENTER',
             'EMB',
         ...

        The result is a dictionary with target code keys and name lists
        as values.  The last name in each list is the one that Skyfield
        uses when printing information about a body.

        """
        d = defaultdict(list)
        for code, name in target_name_pairs:
            if code in self.codes:
                d[code].append(name)
        return dict(d)

    def decode(self, name):
        """Translate a target name into its integer code.

        >>> planets.decode('Venus')
        299

        Raises ``ValueError`` if you supply an unknown name, or
        ``KeyError`` if the target is missing from this kernel.  You can
        supply an integer code if you already have one and just want to
        check whether it is present in this kernel.

        """
        if isinstance(name, int):
            code = name
        else:
            name = name.upper()
            code = _targets.get(name)
            if code is None:
                raise ValueError('unknown SPICE target {0!r}'.format(name))
        if code not in self.codes:
            targets = ', '.join(_format_code_and_name(c) for c in self.codes)
            raise KeyError('kernel {0!r} is missing {1!r} -'
                           ' the targets it supports are: {2}'
                           .format(self.filename, name, targets))
        return code

    def get(self, target):
        # Work in progress: what would it be like if users were not
        # returned Body objects, but general vector objects that they
        # could add and subtract?  Returning raw segments to them is the
        # first step.

        target = self.decode(target)
        segments = self.segments
        segment_dict = dict((segment.target, segment) for segment in segments)
        chain = tuple(_center(target, segment_dict))[::-1]
        if len(chain) == 1:
            return chain[0]
        return Sum(chain[0].center, chain[-1].target, chain)


class VectorFunction(object):
    segments = None

    def __add__(self, other):
        if self.target != other.center:
            if other.target == self.center:
                self, other = other, self
            else:
                raise ValueError()
        return Sum(
            self.center, other.target,
            (self.segments or (self,)) + (other.segments or (other,)),
        )

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        # this should be the only at(), right? and should choose
        # once and for all the right class for the output, Earth etc?
        p, v = self._at(t)
        return build_position(p, v, t, self.center, self.target)


class Segment(VectorFunction):
    __slots__ = ['center', 'target', 'spk_segment']

    def __new__(cls, filename, spk_segment):
        if spk_segment.data_type == 2:
            return object.__new__(ChebyshevPosition)
        if spk_segment.data_type == 3:
            return object.__new__(ChebyshevPositionVelocity)
        raise ValueError('SPK data type {0} not yet supported segment'
                         .format(spk_segment.data_type))

    def __init__(self, filename, spk_segment):
        self.filename = filename
        self.center = spk_segment.center
        self.target = spk_segment.target
        self.spk_segment = spk_segment

    def __str__(self):
        return 'Segment {!r} {}'.format(
            self.filename,
            _format_segment_brief(self),
        )

    def __repr__(self):
        return '<{}>'.format(self)

    def icrf_vector_at(self, t):   # temporary compatibility measure
        return self._at(t)


class ChebyshevPosition(Segment):
    def _at(self, t):
        position, velocity = self.spk_segment.compute_and_differentiate(t.tdb)
        return position / AU_KM, velocity / AU_KM


class ChebyshevPositionVelocity(Segment):
    def _at(self, t):
        pv = self.spk_segment.compute(t.tdb)
        return pv[:3] / AU_KM, pv[3:] * DAY_S / AU_KM


class Sum(VectorFunction):
    def __init__(self, center, target, segments):
        self.center = center
        self.target = target
        self.segments = segments
        self.first = segments[0]
        self.rest = segments[1:]

    def __str__(self):
        segments = self.segments
        lines = '\n'.join('  ' + str(segment) for segment in segments)
        return 'Sum of {} segments:\n{}'.format(len(segments), lines)

    def __repr__(self):
        return '<Sum of {}>'.format(' '.join(repr(s) for s in self.segments))

    def _at(self, t):
        p, v = self.first._at(t)
        for segment in self.rest:
            p2, v2 = segment._at(t)
            p += p2
            v += v2
        return p, v


class Body(object):
    """A target body from a SPICE .bsp kernel file.

    Skyfield programmers usually ask a kernel object to look up and
    return a body object for them, instead of trying to instantiate this
    class directly:

    >>> planets = load('de421.bsp')
    >>> planets['ssb']
    <Body 0 'SOLAR SYSTEM BARYCENTER' from kernel 'de421.bsp'>
    >>> planets[299]
    <Body 299 'VENUS' from kernel 'de421.bsp'>

    """
    def __init__(self, ephemeris, code):
        self.ephemeris = ephemeris
        self.segments = ephemeris.segments
        self.code = code

    def __repr__(self):
        return '<Body {0} from kernel {1!r}>'.format(
            _format_code_and_name(self.code),
            self.ephemeris.filename,
        )

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        """Compute a `Barycentric` position for time `t`.

        The time `t` should be a `Time` object.  The returned position
        will also offer a velocity, if the kernel supports it.

        """
        segments = self.segments
        segment_dict = dict((segment.target, segment) for segment in segments)
        chain = list(_center(self.code, segment_dict))[::-1]
        pos, vel = _tally((), chain, t)
        barycentric = Barycentric(pos, vel, t)
        barycentric.ephemeris = self.ephemeris
        return barycentric

    def geometry_of(self, body):
        """Return a `Geometry` path to another body.

        Given either a `Body` object, or else the name or integer code
        identifying a body in the same ephemeris as this one, compute
        the minimum number of segments necessary to determine their
        relative position and return a `Geometry` object.

        >>> g = earth.geometry_of(moon)
        >>> print(g)
        Geometry from center 399 to target 301 using:
              3 -> 399  EARTH BARYCENTER -> EARTH
              3 -> 301  EARTH BARYCENTER -> MOON

        """
        if not isinstance(body, Body):
            code = self.ephemeris.decode(body)
            body = Body(self.ephemeris, code)
        center_chain, target_chain = _connect(self, body)
        return Geometry(self.code, body.code, center_chain, target_chain)

    def _observe_from_bcrs(self, observer):
        return observe(observer, self)

    def topos(self, latitude=None, longitude=None, latitude_degrees=None,
              longitude_degrees=None, elevation_m=0.0, x=0.0, y=0.0):
        """Return a `Topos` representing a place on Earth.

        See the `Topos` class for a description of the parameters.

        """
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

    def __call__(self, jd):
        """Deprecated alternative to the new at() method."""
        raise DeprecationError("""use method body.at(t), not the call body(t)

If you simply want your old Skyfield script to start working again,
downgrade to Skyfield version 0.4 using a command like:

        pip install skyfield==0.4

Otherwise, you can upgrade your script to modern Skyfield by finding
each place you called a body like a function to generate a position:

        position = body(t)

Instead, Skyfield now offers a method named at(t) to makes the
operation easier to read and more symmetrical with other method calls:

        position = body.at(t)

More documentation can be found at: http://rhodesmill.org/skyfield/""")

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
    """The kernel segments for predicting two bodies' relative position.

    Computing an instantaneous geometry can be faster than computing a
    normal astrometric observation because, instead of referencing both
    bodies back to the Solar System barycenter, the geometry only needs
    the segments that link them.

    For example, the more expensive ``earth.at(t).observe(moon)``
    operation will first compute an Earth position using two kernel
    segments::

        0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
        3 -> 399  EARTH BARYCENTER -> EARTH

    Then it will repeatedly compute the Moon's position as it
    narrows down the light travel time, using two segments::

        0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
        3 -> 301  EARTH BARYCENTER -> MOON

    But a geometry, because it can ignore the real physics required for
    light from one body to reach another, can take a shortcut and avoid
    the Solar System barycenter.  It can also skip the iteration
    required to find the light travel time.

    >>> g = earth.geometry_of(moon)
    >>> print(g)
    Geometry from center 399 to target 301 using:
          3 -> 399  EARTH BARYCENTER -> EARTH
          3 -> 301  EARTH BARYCENTER -> MOON

    Instantaneous geometry positions can be appropriate when plotting a
    diagram of the solar system or of a particular planetary system from
    the point of view of a distant observer.

    """
    def __init__(self, center, target, center_chain, target_chain):
        self.center = center
        self.target = target
        self.center_chain = center_chain
        self.target_chain = target_chain

    def __str__(self):
        segments = self.center_chain + self.target_chain
        lines = '\n'.join(_format_segment(s) for s in segments)
        return 'Geometry from center {0} to target {1} using:\n{2}'.format(
            self.center, self.target, lines)

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        """Compute instantaneous position and velocity between two bodies.

        The argument ``t`` should be a `Time`, and the return value will
        be the relative position of the the target body relative to the
        center body.  If the center is 0 "Solar System Barycenter" then
        the result will be `Barycentric`.  Otherwise, it will be a plain
        `ICRF` position.

        """
        pos, vel = _tally(self.center_chain, self.target_chain, t)
        cls = Barycentric if self.center == 0 else ICRF
        return cls(pos, vel, t)


def _format_code_and_name(code):
    name = _names.get(code, None)
    if name is None:
        return str(code)
    return '{0} {1!r}'.format(code, name)

def _format_segment(segment):
    cname = _names.get(segment.center, 'unknown')
    tname = _names.get(segment.target, 'unknown')
    return '    {0:3} -> {1:<3}  {2} -> {3}'.format(
        segment.center, segment.target, cname, tname)

def _format_segment_brief(segment):
    cname = _names.get(segment.center)
    tname = _names.get(segment.target)
    return '{}{}{} -> {}{}{}'.format(
        segment.center,
        ' ' if cname else '',
        cname,
        segment.target,
        ' ' if tname else '',
        tname,
    )

def _tally(minus_chain, plus_chain, t):
    position = velocity = 0.0
    for segment in minus_chain:
        p, v = segment.icrf_vector_at(t)
        position -= p
        velocity -= v
    for segment in plus_chain:
        p, v = segment.icrf_vector_at(t)
        position += p
        velocity += v
    return position, velocity
