"""Open a BPC file, read its angles, and produce rotation matrices."""

import re
from .functions import rot_x, rot_z
from .units import Angle

class PlanetaryConstants(object):
    """Planetary constants kernel."""

    def __init__(self):
        self.assignments = {}
        self.segments = {}

    def load(self, something):
        """Yeah."""

    def load_text_pck(self, lines):
        self.assignments.update(parse_text_pck(lines))

    def load_binary_pck(self):
        pass

class PlanetaryConstantsSegment(object):
    def rotation_at(self, t):
        ra, dec, w = eph.segments[0].compute(tdb, 0.0, False)
        return einsum('ij...,jk...,kl...->il...',
                      rot_z(-w), rot_x(-dec), rot_z(-ra))

def parse_text_pck(lines):
    """Yield ``(name, value)`` tuples parsed from a PCK text kernel."""
    tokens = iter(_parse_text_pck_tokens(lines))
    for token in tokens:
        name = token
        equals = next(tokens)
        if equals != '=':
            raise ValueError('was expecting an equals sign after %r' % name)
        value = next(tokens)
        if value == '(':
            value = []
            for token2 in tokens:
                if token2 == ')':
                    break
                value.append(_evaluate(token2))
        else:
            value = _evaluate(value)
        yield name, value

def _evaluate(token):
    """Return a string, integer, or float parsed from a PCK text kernel."""
    if token[0] == "'":
        return token[1:-1]
    if token.isdigit():
        return int(token)
    if token.startswith('@'):
        raise NotImplemented('TODO: need parser for dates, like @01-MAY-1991/16:25')
    return float(token)

_token_re = re.compile(r"'[^']*'|[^', ]+")

def _parse_text_pck_tokens(lines):
    """Yield all the tokens inside the data segments of a PCK text file."""
    lines = iter(lines)
    for line in lines:
        line = line.strip()
        if line != r'\begindata':
            continue
        for line in lines:
            line = line.strip()
            if line == r'\begintext':
                break
            for token in _token_re.findall(line):
                yield token
