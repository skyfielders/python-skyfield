"""An interface between JPL ephemerides and Skyfield."""

from collections import namedtuple
from numpy import max, min

from jplephem.spk import SPK
from jplephem.names import target_names as _names

from .constants import AU_KM, C_AUDAY, DAY_S
from .errors import DeprecationError, raise_error_for_deprecated_time_arguments
from .functions import length_of
from .positionlib import Astrometric, Barycentric, ICRF
from .timelib import calendar_date

Segment = namedtuple('Segment', 'center target compute')
_targets = dict((name, target) for (target, name) in _names.items())


class SpiceKernel(object):
    """Ephemeris file from NASA's SPICE system

    You can download a .bsp file yourself and use this class to open it,
    or use the Skyfield `load()` function to find a popular ephemeris.
    Once loaded, you can print a SPICE kernel to the screen to see the
    Solar System bodies for which it can generate positions:

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

    You can look up objects by name or by their integer code:

    >>> planets['earth']
    <Body EARTH from kernel 'de421.bsp'>
    >>> planets[499]
    <Body MARS from kernel 'de421.bsp'>

    """
    def __init__(self, filename):
        self.filename = filename
        self.spk = SPK.open(filename)
        self.segments = [Segment(s.center, s.target, _build_compute(s))
                         for s in self.spk.segments]
        self.codes = set(s.center for s in self.segments).union(
                         s.target for s in self.segments)

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
            cname = _names.get(s.center, 'unknown')
            tname = _names.get(s.target, 'unknown')
            lines.append('    {0:3} -> {1:<3}  {2} -> {3}'
                         .format(s.center, s.target, cname, tname))
        return '\n'.join(lines)

    def __getitem__(self, name):
        """Return a `Body` for a SPICE target name or integer"""
        code = self.decode(name)
        return Body(self, code)

    def decode(self, name):
        """Given a SPICE target name, return its integer"""
        if isinstance(name, int):
            return name
        name = name.upper()
        code = _targets.get(name)
        if code is None:
            raise KeyError('unknown SPICE target name {0!r}'.format(name))
        if code not in self.codes:
            names = ', '.join(_names[c] for c in self.codes)
            raise KeyError('kernel {0!r} is missing {1!r} -'
                           ' the targets it supports are: {2}'
                           .format(self.filename, name, names))
        return code


def _build_compute(segment):
    """Build a Skyfield `compute` callback for the SPK `segment`."""

    if segment.data_type == 2:
        def compute(t):
            position, velocity = segment.compute_and_differentiate(t.tdb)
            return position / AU_KM, velocity / AU_KM

    elif segment.data_type == 3:
        def compute(t):
            six = segment.compute(t.tdb)
            return six[:3] / AU_KM, six[3:] * DAY_S / AU_KM

    else:
        raise ValueError('SPK data type {0} not yet supported segment'
                         .format(segment.data_type))
    return compute


class Body(object):
    def __init__(self, ephemeris, code):
        self.ephemeris = ephemeris
        self.segments = ephemeris.segments
        self.code = code

    def __repr__(self):
        return '<Body {0} from kernel {1!r}>'.format(
            _names.get(self.code, self.code), self.ephemeris.filename)

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        """Compute the `Barycentric` position at time `t`"""
        segments = self.segments
        segment_dict = dict((segment.target, segment) for segment in segments)
        chain = list(_center(self.code, segment_dict))[::-1]
        pos, vel = _tally((), chain, t)
        barycentric = Barycentric(pos, vel, t)
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
        raise ValueError('observe_from() light-travel time'
                         ' failed to converge')
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
    def __init__(self, center, target, center_chain, target_chain):
        self.center = center
        self.target = target
        self.center_chain = center_chain
        self.target_chain = target_chain

    def __str__(self):
        return 'Geometry\n{0}'.format('\n'.join(
            ' {0}'.format(c)
            for c in self.center_chain + self.target_chain))

    @raise_error_for_deprecated_time_arguments
    def at(self, t):
        """Return the geometric Cartesian position and velocity."""
        pos, vel = _tally(self.center_chain, self.target_chain, t)
        cls = Barycentric if self.center == 0 else ICRF
        return cls(pos, vel, t)


def _tally(minus_chain, plus_chain, t):
    position = velocity = 0.0
    for segment in minus_chain:
        p, v = segment.compute(t)
        position -= p
        velocity -= v
    for segment in plus_chain:
        p, v = segment.compute(t)
        position += p
        velocity += v
    return position, velocity
