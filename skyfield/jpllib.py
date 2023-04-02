"""An interface between JPL ephemerides and Skyfield."""

import os
from collections import defaultdict

from jplephem.exceptions import OutOfRangeError
from jplephem.spk import SPK
from jplephem.names import target_name_pairs

from .constants import AU_KM, DAY_S
from .errors import EphemerisRangeError
from .timelib import compute_calendar_date
from .vectorlib import VectorFunction, VectorSum, _jpl_code_name_dict

_jpl_name_code_dict = dict(
    (name, target) for (target, name) in target_name_pairs
)

class SpiceKernel(object):
    """Ephemeris file in NASA .bsp format.

    A "Spacecraft and Planet Kernel" (SPK) file from NASA provides
    |xyz| coordinates for bodies in the Solar System like the Sun,
    planets, moons, and spacecraft.

    You can download a .bsp file yourself and use this class to open it,
    or use the Skyfield ``load()`` function to automatically download a
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

    To retrieve the one or more vectors necessary to compute the
    position of a body relative to the Solar System barycenter, look up
    the body by its name or official SPICE identifying integer:

    >>> planets['earth']
    <VectorSum of 2 vectors:
     'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 3 EARTH BARYCENTER
     'de421.bsp' segment 3 EARTH BARYCENTER -> 399 EARTH>
    >>> planets[499]
    <VectorSum of 2 vectors:
     'de421.bsp' segment 0 SOLAR SYSTEM BARYCENTER -> 4 MARS BARYCENTER
     'de421.bsp' segment 4 MARS BARYCENTER -> 499 MARS>

    The result will be a :class:`~skyfield.vectorlib.VectorFunction`
    instance that you can ask for a position at a given input time.

    """
    def __init__(self, path):
        self.path = path
        self.filename = os.path.basename(path)
        self.spk = SPK.open(path)
        self.segments = [SPICESegment(self, s) for s in self.spk.segments]
        self.codes = set(s.center for s in self.segments).union(
                         s.target for s in self.segments)

    def __repr__(self):
        return '<{0} {1!r}>'.format(type(self).__name__, self.path)

    def __str__(self):
        segments = self.spk.segments
        lines = ['SPICE kernel file {0!r} has {1} segments'
                 .format(self.filename, len(segments))]
        format_date = '{0}-{1:02}-{2:02}'.format
        start = end = None
        for s in segments:
            if start != s.start_jd or end != s.end_jd:
                start, end = s.start_jd, s.end_jd
                starts = format_date(*compute_calendar_date(int(start)))
                ends = format_date(*compute_calendar_date(int(end)))
                lines.append('  JD {0:.2f} - JD {1:.2f}  ({2} through {3})'
                             .format(start, end, starts, ends))
            lines.append(_format_segment(s))
        return '\n'.join(lines)

    def close(self):
        """Close this ephemeris file."""
        self.spk.close()

        # In practice, users are not confident the file is really closed
        # unless the metadata also disappears.
        del self.segments[:]
        self.codes.clear()

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
            code = _jpl_name_code_dict.get(name)
            if code is None:
                raise ValueError('unknown SPICE target {0!r}'.format(name))
        if code not in self.codes:
            targets = ', '.join(_format_code_and_name(c) for c in self.codes)
            raise KeyError('kernel {0!r} is missing {1!r} -'
                           ' the targets it supports are: {2}'
                           .format(self.filename, name, targets))
        return code

    def __getitem__(self, target):
        """Return a vector function for computing the location of `target`."""
        target = self.decode(target)
        segments = self.segments
        segment_dict = dict((segment.target, segment) for segment in segments)
        chain = tuple(_center(target, segment_dict))
        if len(chain) == 1:
            return chain[0]
        chain = chain[::-1]
        center = chain[0].center
        target = chain[-1].target
        return VectorSum(center, target, chain)

    def __contains__(self, name_or_code):
        if isinstance(name_or_code, int):
            code = name_or_code
        else:
            code = _jpl_name_code_dict.get(name_or_code.upper())
        return code in self.codes

class SPICESegment(VectorFunction):

    def __new__(cls, ephemeris, spk_segment):
        if spk_segment.data_type == 2:
            return object.__new__(ChebyshevPosition)
        if spk_segment.data_type == 3:
            return object.__new__(ChebyshevPositionVelocity)
        raise ValueError('SPK data type {0} not yet supported'
                         .format(spk_segment.data_type))

    def __init__(self, ephemeris, spk_segment):
        self.ephemeris = ephemeris
        self.center = spk_segment.center
        self.target = spk_segment.target
        self.spk_segment = spk_segment

    @property
    def vector_name(self):
        return '{0!r} segment'.format(self.ephemeris.path)

    def time_range(self, ts):
        s = self.spk_segment
        return ts.tdb_jd(s.start_jd), ts.tdb_jd(s.end_jd)

class ChebyshevPosition(SPICESegment):
    def _at(self, t):
        segment = self.spk_segment
        try:
            position, velocity = segment.compute_and_differentiate(
                t.whole, t.tdb_fraction)
        except OutOfRangeError as e:
            start_time, end_time = self.time_range(t.ts)
            s = '%04d-%02d-%02d' % start_time.tdb_calendar()[:3]
            t = '%04d-%02d-%02d' % end_time.tdb_calendar()[:3]
            text = 'ephemeris segment only covers dates %s through %s' % (s, t)
            mask = e.out_of_range_times
            segment = self.spk_segment
            e = EphemerisRangeError(text, start_time, end_time, mask, segment)
            e.__cause__ = None  # avoid exception chaining in Python 3
            raise e

        return position / AU_KM, velocity / AU_KM, None, None


class ChebyshevPositionVelocity(SPICESegment):
    def _at(self, t):
        pv = self.spk_segment.compute(t.whole, t.tdb_fraction)
        return pv[:3] / AU_KM, pv[3:] * DAY_S / AU_KM, None, None


def _center(code, segment_dict):
    """Starting with `code`, follow segments from target to center."""
    while code in segment_dict:
        segment = segment_dict[code]
        yield segment
        code = segment.center

def _format_code_and_name(code):
    name = _jpl_code_name_dict.get(code, None)
    if name is None:
        return str(code)
    return '{0} {1}'.format(code, name)

def _format_segment(segment):
    cname = _jpl_code_name_dict.get(segment.center, 'unknown')
    tname = _jpl_code_name_dict.get(segment.target, 'unknown')
    return '    {0:3} -> {1:<3}  {2} -> {3}'.format(
        segment.center, segment.target, cname, tname)
