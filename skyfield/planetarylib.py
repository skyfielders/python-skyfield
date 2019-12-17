"""Open a BPC file, read its angles, and produce rotation matrices."""

import re
from numpy import array, einsum, nan
from jplephem.pck import DAF, PCK
from .constants import ASEC2RAD, AU_KM
from .functions import rot_x, rot_y, rot_z
from .units import Angle, Distance
from .vectorlib import VectorFunction

_TEXT_MAGIC_NUMBERS = b'KPL/FK', b'KPL/PCK'
_NAN3 = array((nan, nan, nan))

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
        """Do .assignments[key] but with a pretty exception on failure."""
        try:
            return self.assignments[key]
        except KeyError:
            e = ValueError(_missing_name_message.format(key))
            e.__cause__ = None
            raise e

    def build_frame_named(self, name):
        integer = self._get_assignment('FRAME_{0}'.format(name))
        return self.build_frame(integer)

    def build_frame(self, integer, _segment=None):
        center = self._get_assignment('FRAME_{0}_CENTER'.format(integer))
        spec = self.assignments.get('TKFRAME_{0}_SPEC'.format(integer))
        if spec is None:
            matrix = None
        else:
            if spec == 'ANGLES':
                angles = self.assignments['TKFRAME_{0}_ANGLES'.format(integer)]
                axes = self.assignments['TKFRAME_{0}_AXES'.format(integer)]
                units = self.assignments['TKFRAME_{0}_UNITS'.format(integer)]
                scale = _unit_scales[units]
                matrix = 1,0,0, 0,1,0, 0,0,1
                matrix = array(matrix)
                matrix.shape = 3, 3
                # TODO: test is not yet sensitive enough to know if we
                # did these rotations in the right order.
                # (Can we reverse order without loss of correctness?)
                for angle, axis in list(zip(angles, axes)):
                    rot = _rotations[axis]
                    #matrix = rot(angle * scale).dot(matrix)
                    matrix = matrix.dot(rot(angle * scale))
            elif spec == 'MATRIX':
                matrix = self.assignments['TKFRAME_{0}_MATRIX'.format(integer)]
                matrix = array(matrix)
                matrix.shape = 3, 3
            else:
                raise NotImplemented('spec %r not yet implemented' % spec)
            relative = self.assignments['TKFRAME_{0}_RELATIVE'.format(integer)]
            integer = self.assignments['FRAME_{0}'.format(relative)]

        if _segment is None:
            segment = self._segment_map.get(integer)
        else:
            segment = _segment

        if segment is None:
            raise LookupError('you have not yet loaded a binary PCK file that'
                              ' has a segment for frame {0}'.format(integer))
        assert segment.frame == 1  # base frame should be ITRF/J2000
        return Frame(center, segment, matrix)

    def build_latlon_degrees(self, frame, latitude_degrees, longitude_degrees):
        lat = Angle.from_degrees(latitude_degrees)
        lon = Angle.from_degrees(longitude_degrees)
        radii = self._get_assignment('BODY{0}_RADII'.format(frame.center))
        if not radii[0] == radii[1] == radii[2]:
            raise ValueError('only spherical bodies are supported,'
                             ' but the radii of this body are: %s' % radii)
        distance = Distance(au=radii[0] / AU_KM)
        return PlanetTopos.from_latlon_distance(frame, lat, lon, distance)

_rotations = None, rot_x, rot_y, rot_z
_unit_scales = {'ARCSECONDS': ASEC2RAD}
_missing_name_message = """unknown planetary constant {0!r}

You should either use this object's `.read_text()` method to load an
additional "*.tf" PCK text file that defines the missing name, or
manually provide a value by adding the name and value to the this
object's `.assignments` dictionary."""

class Frame(object):
    """Planetary constants frame, for building rotation matrices."""

    def __init__(self, center, segment, matrix):
        self.center = center
        self._segment = segment
        self._matrix = matrix

    def rotation_at(self, t):
        ra, dec, w = self._segment.compute(t.tdb, 0.0, False)
        R = rot_z(-w).dot(rot_x(-dec).dot(rot_z(-ra)))
        if self._matrix is not None:
            R = self._matrix.dot(R)
        return R

class PlanetTopos(VectorFunction):
    """Location that rotates with the surface of another Solar System body.

    The location can either be on the surface of the body, or in some
    other fixed position that rotates with the body's surface.

    """
    def __init__(self, center, frame, position_au):
        # TODO: always take center from frame
        self.center = center
        self.target = object()  # TODO: make more interesting
        self.center_name = None  # TODO: deprecate and remove
        self.target_name = None
        self._frame = frame
        self._position_au = position_au

    @classmethod
    def from_latlon_distance(cls, frame, latitude, longitude, distance):
        r = array((-distance.au, 0.0, 0.0))
        r = rot_z(longitude.radians).dot(rot_y(-latitude.radians).dot(r))

        self = cls(frame.center, frame, r)
        self.latitude = latitude
        self.longitude = longitude
        return self

    def _at(self, t):
        r = self._frame.rotation_at(t).T.dot(self._position_au)
        v = _NAN3.copy()
        # TODO: altaz
        return r, v, r, None

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

_token_re = re.compile(b"'[^']*'|[^', ]+")

def _parse_text_pck_tokens(lines):
    """Yield all the tokens inside the data segments of a PCK text file."""
    lines = iter(lines)
    for line in lines:
        if b'\\begindata' not in line:
            continue  # save cost of strip() on most lines
        line = line.strip()
        if line != b'\\begindata':
            continue
        for line in lines:
            line = line.strip()
            if line == b'\\begintext':
                break
            for token in _token_re.findall(line):
                yield token
