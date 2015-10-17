"""n interface between JPL ephemerides and Skyfield."""

from jplephem.spk import SPK
from jplephem.names import target_names as _names

from .constants import AU_KM, DAY_S
from .ephemerislib import Body, Segment

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
