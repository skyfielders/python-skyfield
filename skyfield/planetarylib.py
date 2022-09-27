# -*- coding: utf-8 -*-
"""Open a BPC file, read its angles, and produce rotation matrices."""

from numpy import array, cos, nan, sin, deg2rad
from jplephem.pck import DAF, PCK
from .constants import ASEC2RAD, AU_KM, DAY_S, tau, T0, D_JC
from .data import text_pck
from .functions import _T, mxv, mxm, mxmxm, rot_x, rot_y, rot_z
from .units import Angle, Distance
from .vectorlib import VectorFunction

_TEXT_MAGIC_NUMBERS = b'KPL/FK', b'KPL/PCK'
_NAN3 = array((nan, nan, nan))
_halftau = tau / 2.0
_quartertau = tau / 4.0

class PlanetaryConstants(object):
    """Planetary constants manager.

    You can use this class to build working models of Solar System
    bodies by loading text planetary constant files and binary
    orientation kernels.  For a full description of how to use this, see
    :doc:`planetary`.

    """
    def __init__(self):
        self.variables = {}
        self._binary_files = []
        self._segment_map = {}

    @property
    def assignments(self):  # deprecated original name for the variables dict
        return self.variables

    def read_text(self, file):
        """Read frame variables from a KPL/FK file.

        Appropriate files will typically have the extension ``.tf`` or
        ``.tpc`` and will define a series of names and values that will
        be loaded into this object's ``.variables`` dictionary.

        >>> from skyfield.api import load
        >>> pc = PlanetaryConstants()
        >>> pc.read_text(load('moon_080317.tf'))
        >>> pc.variables['FRAME_31006_NAME']
        'MOON_PA_DE421'

        """
        file.seek(0)
        try:
            if not file.read(7).startswith(_TEXT_MAGIC_NUMBERS):
                raise ValueError('file must start with one of the patterns:'
                                 ' {0}'.format(_TEXT_MAGIC_NUMBERS))
            text_pck.load(file, self.variables)
        finally:
            file.close()

    def read_binary(self, file):
        """Read binary segments descriptions from a DAF/PCK file.

        Binary segments live in ``.bpc`` files and predict how a body
        like a planet or moon will be oriented on a given date.

        """
        file.seek(0)
        if file.read(7) != b'DAF/PCK':
            raise ValueError('file must start with the bytes "DAF/PCK"')
        pck = PCK(DAF(file))
        self._binary_files.append(pck)
        for segment in pck.segments:
            self._segment_map[segment.body] = segment

    def _get_assignment(self, key):
        """Do .variables[key] but with a pretty exception on failure."""
        try:
            return self.variables[key]
        except KeyError:
            e = ValueError(_missing_name_message.format(key))
            e.__cause__ = None
            raise e

    def build_frame_named(self, name):
        """Given a frame name, return a :class:`Frame` object."""
        integer = self._get_assignment('FRAME_{0}'.format(name))
        return self.build_frame(integer)

    def build_frame(self, integer, _segment=None):
        """Given a frame integer code, return a :class:`Frame` object."""
        center = self._get_assignment('FRAME_{0}_CENTER'.format(integer))
        spec = self.variables.get('TKFRAME_{0}_SPEC'.format(integer))
        if spec is None:
            matrix = None
        else:
            if spec == 'ANGLES':
                angles = self.variables['TKFRAME_{0}_ANGLES'.format(integer)]
                axes = self.variables['TKFRAME_{0}_AXES'.format(integer)]
                units = self.variables['TKFRAME_{0}_UNITS'.format(integer)]
                scale = _unit_scales[units]
                matrix = 1,0,0, 0,1,0, 0,0,1
                matrix = array(matrix)
                matrix.shape = 3, 3
                for angle, axis in list(zip(angles, axes)):
                    rot = _rotations[axis]
                    matrix = mxm(rot(angle * scale), matrix)
            elif spec == 'MATRIX':
                matrix = self.variables['TKFRAME_{0}_MATRIX'.format(integer)]
                matrix = array(matrix)
                matrix.shape = 3, 3
            else:
                raise NotImplementedError('spec %r not yet implemented' % spec)
            relative = self.variables['TKFRAME_{0}_RELATIVE'.format(integer)]
            integer = self.variables['FRAME_{0}'.format(relative)]

        if _segment is None:
            segment = self._segment_map.get(integer)
        else:
            segment = _segment

        if segment is None:
            raise LookupError('you have not yet loaded a binary PCK file that'
                              ' has a segment for frame {0}'.format(integer))
        assert segment.frame == 1  # base frame should be ITRF/J2000
        return Frame(center, segment, matrix)

    def build_latlon_degrees(self, frame, latitude_degrees, longitude_degrees,
                             elevation_m=0.0):
        """Build an object representing a location on a body's surface."""
        lat = Angle.from_degrees(latitude_degrees)
        lon = Angle.from_degrees(longitude_degrees)
        radii = self._get_assignment('BODY{0}_RADII'.format(frame.center))
        if not radii[0] == radii[1] == radii[2]:
            raise ValueError('only spherical bodies are supported,'
                             ' but the radii of this body are: %s' % radii)
        au = (radii[0] + elevation_m * 1e-3) / AU_KM
        distance = Distance(au)
        return PlanetTopos.from_latlon_distance(frame, lat, lon, distance)

_rotations = None, rot_x, rot_y, rot_z
_unit_scales = {'ARCSECONDS': ASEC2RAD}
_missing_name_message = """unknown planetary constant {0!r}

You should either use this object's `.read_text()` method to load an
additional "*.tf" PCK text file that defines the missing name, or
manually provide a value by adding the name and value to the this
object's `.variables` dictionary."""

class Frame(object):
    """Planetary constants frame, for building rotation matrices."""

    def __init__(self, center, segment, matrix):
        self.center = center
        self._segment = segment
        self._matrix = matrix

    def rotation_at(self, t):
        """Return the rotation matrix for this frame at time ``t``."""
        ra, dec, w = self._segment.compute(t.tdb, 0.0, False)
        R = mxm(rot_z(-w), mxm(rot_x(-dec), rot_z(-ra)))
        if self._matrix is not None:
            R = mxm(self._matrix, R)
        return R

    def rotation_and_rate_at(self, t):
        """Return rotation and rate matrices for this frame at time ``t``.

        The rate matrix returned is in units of angular motion per day.

        """
        components, rates = self._segment.compute(t.whole, t.tdb_fraction, True)
        ra, dec, w = components
        radot, decdot, wdot = rates

        R = mxm(rot_z(-w), mxm(rot_x(-dec), rot_z(-ra)))

        zero = w * 0.0
        ca = cos(w)
        sa = sin(w)
        u = cos(dec)
        v = -sin(dec)

        domega0 = wdot + u * radot
        domega1 = ca * decdot - sa * v * radot
        domega2 = sa * decdot + ca * v * radot

        drdtrt = array((
            (zero, domega0, domega2),
            (-domega0, zero, domega1),
            (-domega2, -domega1, zero),
        ))

        dRdt = mxm(drdtrt, R)

        if self._matrix is not None:
            R = mxm(self._matrix, R)
            dRdt = mxm(self._matrix, dRdt)

        return R, dRdt * DAY_S

class PlanetTopos(VectorFunction):
    """Location that rotates with the surface of another Solar System body.

    The location can either be on the surface of the body, or in some
    other fixed position that rotates with the body's surface.

    """
    def __init__(self, frame, position_au):
        self.center = frame.center
        self._frame = frame
        self._position_au = position_au

    @classmethod
    def from_latlon_distance(cls, frame, latitude, longitude, distance):
        r = array((distance.au, 0.0, 0.0))
        r = mxv(rot_z(longitude.radians), mxv(rot_y(-latitude.radians), r))

        self = cls(frame, r)
        self.latitude = latitude
        self.longitude = longitude
        return self

    @property
    def target(self):
        # When used as a vector function, this planetary geographic
        # location computes positions from the planet's center to
        # itself.  (This is a property, rather than an attribute, to
        # avoid a circular reference that delays garbage collection.)
        return self

    def _at(self, t):
        # Since `_position_au` has zero velocity in this reference
        # frame, velocity includes a `dRdt` term but not an `R` term.
        R, dRdt = self._frame.rotation_and_rate_at(t)
        r = mxv(_T(R), self._position_au)
        v = mxv(_T(dRdt), self._position_au)
        return r, v, None, None

    def rotation_at(self, t):
        """Compute the altazimuth rotation matrix for this location’s sky."""
        R = mxmxm(
            # TODO: Figure out how to produce this rotation directly
            # from _position_au, to support situations where we were not
            # given a latitude and longitude.  If that is not feasible,
            # then at least cache the product of these first two matrices.
            rot_y(_quartertau - self.latitude.radians),
            rot_z(_halftau - self.longitude.radians),
            self._frame.rotation_at(t),
        )
        # TODO:
        # Can clockwise be turned into counterclockwise through any
        # possible rotation?  For now, flip the sign of y so that
        # azimuth reads north-east rather than the other direction.
        R[1] *= -1
        return R

class PlanetaryOrientation():
  """Planet/moon orientation and shape using models from .tpc file
  
  Based on the models presented in pck00010.tpc, which is based on
     
  Archinal, B.A., A’Hearn, M.F., Bowell, E. et al. Report of the IAU
  Working Group on Cartographic Coordinates and Rotational Elements: 2009.
  Celest Mech Dyn Astr 109, 101–135 (2011).
  https://doi.org/10.1007/s10569-010-9320-4

  >>> from skyfield.api import load, PlanetaryConstants, PlanetaryOrientation
  >>> pc = PlanetaryConstants()
  >>> pc.read_text(load('pck00010.tpc'))
  >>> jup_orientation = PlanetaryOrientation(pc, '599')
  >>> t = load.timescale().utc(2020,1,1)
  >>> jup_orientation.at(t).orientation # [RA, DEC, prime_meridian]
  [268.05773343067165, 64.49711711016498, 6359115.859072568]

  """
  def __init__(self, pc, naifid):
    self.a0 = None
    self.d0 = None
    self.W = None
    self.orientation = [self.a0, self.d0, self.W]

    pck = pc.variables

    # make sure that naifid is a string number
    naifid = str(naifid)
    if naifid.upper().startswith('BODY'):
      naifid = naifid[4:]

    self.naifid = naifid

    # not all variables are available for every body
    # set to None if they aren't
    self.pole_ra = self._getPCKvar(pck, naifid, 'POLE_RA')
    self.pole_dec = self._getPCKvar(pck, naifid, 'POLE_DEC')
    self.pm = self._getPCKvar(pck, naifid, 'PM')
    self.long_axis = self._getPCKvar(pck, naifid, 'LONG_AXIS')
    self.nut_prec_ra = self._getPCKvar(pck, naifid, 'NUT_PREC_RA')
    self.nut_prec_dec = self._getPCKvar(pck, naifid, 'NUT_PREC_DEC')
    self.nut_prec_pm = self._getPCKvar(pck, naifid, 'NUT_PREC_PM')
    # NUT_PREC_ANGLES are listed under main planet barycenter number
    self.nut_prec_angles = self._getPCKvar(pck, naifid[0], 'NUT_PREC_ANGLES')
    self.radii = self._getPCKvar(pck, naifid, 'RADII')

    # check if POLE and PM values (orientation) are available
    # many moons don't have them
    if (self.pole_ra is not None and
        self.pole_dec is not None and
        self.pm is not None):
      self.orientation_model = True
    else:
      self.orientation_model = False

    # check if _NUT_PREC_* values (nutation) are available
    # not all bodies have it, most notably Earth
    if (self.nut_prec_ra is not None and
        self.nut_prec_dec is not None and
        self.nut_prec_pm is not None and
        self.nut_prec_angles is not None):
      self.nutation_model = True
    else:
      self.nutation_model = False

    # check if RADII values (shape) are available
    if (self.radii is not None):
      self.shape_model = True
    else:
      self.shape_model = False

    # check that at least one model is available
    if (not self.orientation_model and
        not self.nutation_model and
        not self.shape_model):
      e = 'No models found for BODY' + self.naifid + \
          '. Verify that NAIF ID is correct'
      raise KeyError(e)

    # parse NUT_PREC_ANGLES into list of tuples for easier use later
    if self.nutation_model:
      self.nut_prec_angles_pairs = list(zip( self.nut_prec_angles[::2],
                                             self.nut_prec_angles[1::2] ))

  def _getPCKvar(self, pckvars, naifid, var):
    """Do .variables[key] but return None instead of throwing an exception."""
    try:
      key = 'BODY'+naifid+'_'+var # e.g. BODY399_POLE_RA
      v = pckvars[key]
    except KeyError:
      v = None
    return v

  def at(self, t):
    self.t = t
    self._calcOrientation()
    return self

  def _calcOrientation(self, t=None):
    if not self.orientation_model:
      e = 'BODY' + self.naifid + ' does not have an orientation model'
      raise AttributeError(e)

    if t is None:
      t = self.t
    else:
      self.t = t
    
    d = t.tdb - T0 # Interval in days from the standard epoch
    T = d/D_JC # Interval in Julian centuries from the standard epoch

    nut_prec_ra_sum = 0
    nut_prec_dec_sum = 0
    nut_prec_pm_sum = 0

    # only do the nutation calculations if the necessary values are available
    if self.nutation_model:
      # Is keying to NUT_PREC_RA really the correct thing to do?
      # This keeps everything in a single loop and there doesn't seem to be
      # any cases NUT_PREC variables for the same body having different numbers
      # of values.
      for i in range(len(self.nut_prec_ra)):
        a1, a2 = self.nut_prec_angles_pairs[i]
        npa = a1 + a2*T
        npa_rad = deg2rad(npa)

        nut_prec_ra_sum += self.nut_prec_ra[i]*sin(npa_rad)
        nut_prec_dec_sum += self.nut_prec_dec[i]*cos(npa_rad)
        nut_prec_pm_sum += self.nut_prec_pm[i]*sin(npa_rad)

    a0 = self.pole_ra[0] + self.pole_ra[1]*T + nut_prec_ra_sum
    d0 = self.pole_dec[0] + self.pole_dec[1]*T + nut_prec_dec_sum
    W = self.pm[0] + self.pm[1]*d + nut_prec_pm_sum

    self.orientation = [a0, d0, W]
    self.a0, self.d0, self.W = self.orientation

    return self.orientation
