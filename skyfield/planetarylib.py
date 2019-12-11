"""Open a BPC file, read its angles, and produce rotation matrices."""

import re
from numpy import einsum
from jplephem.pck import DAF, PCK
from .functions import rot_x, rot_z
from .units import Angle

_TEXT_MAGIC_NUMBERS = b'KPL/FK', b'KPL/PCK'

class PlanetaryConstants(object):
    """Planetary constants kernel."""

    def __init__(self):
        self.assignments = {}
        self._binary_files = []
        self._segment_map = {}

    def read_text(self, file):
        """Read frame assignments from a KPL/FK file."""
        file.seek(0)
        try:
            if not file.read(7).startswith(_TEXT_MAGIC_NUMBERS):
                raise ValueError('file must start with one of the patterns:'
                                 ' {0}'.format(_TEXT_MAGIC_NUMBERS))
            file.seek(0)
            self.assignments.update(parse_text_pck(file))
        finally:
            file.close()

    def read_binary(self, file):
        """Read binary segments descriptions from a DAF/PCK file."""
        file.seek(0)
        if file.read(7) != b'DAF/PCK':
            raise ValueError('file must start with the bytes "DAF/PCK"')
        pck = PCK(DAF(file))
        self._binary_files.append(pck)
        for segment in pck.segments:
            self._segment_map[segment.body] = segment

    def _get_assignment(self, key):
        try:
            return self.assignments[key]
        except KeyError:
            e = ValueError(_missing_name_message.format(key))
            e.__cause__ = None
            raise e

    def build_frame_named(self, name):
        integer = self._get_assignment('FRAME_{0}'.format(name))
        return self.build_frame(integer)

    def build_frame(self, integer):
        center = self._get_assignment('FRAME_{0}_CENTER'.format(integer))
        segment = self._segment_map[integer]
        assert segment.frame == 1  # base frame should be ITRF/J2000
        return Frame(segment)

_missing_name_message = """unknown planetary constant {0!r}

You should either use this object's `.read_text()` method to load an
additional "*.tf" PCK text file that defines the missing constant, or
manually provide a value by adding the key and value to the this
object's `.assignments` dictionary."""

class Frame(object):
    """Planetary constants frame, for building rotation matrices."""

    def __init__(self, segment):
        self._segment = segment

    def rotation_at(self, t):
        ra, dec, w = self._segment.compute(t.tdb, 0.0, False)
        return einsum('ij...,jk...,kl...->il...',
                      rot_z(-w), rot_x(-dec), rot_z(-ra))

class FramePosition(object):

    def __init__(self, frame):
        self._frame = frame

    def at(self, t):
        R = self._frame.at(t)

def parse_text_pck(lines):
    """Yield ``(name, value)`` tuples parsed from a PCK text kernel."""
    tokens = iter(_parse_text_pck_tokens(lines))
    for token in tokens:
        name = token.decode('ascii')
        equals = next(tokens)
        if equals != b'=':
            raise ValueError('was expecting an equals sign after %r' % name)
        value = next(tokens)
        if value == b'(':
            value = []
            for token2 in tokens:
                if token2 == b')':
                    break
                value.append(_evaluate(token2))
        else:
            value = _evaluate(value)
        yield name, value

def _evaluate(token):
    """Return a string, integer, or float parsed from a PCK text kernel."""
    if token[0:1].startswith(b"'"):
        return token[1:-1].decode('ascii')
    if token.isdigit():
        return int(token)
    if token.startswith(b'@'):
        raise NotImplemented('TODO: need parser for dates, like @01-MAY-1991/16:25')
    token = token.replace(b'D', b'E')  # for numbers like -1.4D-12
    return float(token)

_token_re = re.compile(rb"'[^']*'|[^', ]+")

def _parse_text_pck_tokens(lines):
    """Yield all the tokens inside the data segments of a PCK text file."""
    lines = iter(lines)
    for line in lines:
        line = line.strip()
        if line != rb'\begindata':
            continue
        for line in lines:
            line = line.strip()
            if line == rb'\begintext':
                break
            for token in _token_re.findall(line):
                yield token
