"""Simple distance, velocity, and angle support for Skyfield.

"""
from __future__ import print_function

import numpy as np
import sys
from numpy import abs, array, copysign, isnan
from .constants import AU_KM, AU_M, C, DAY_S, tau
from .descriptorlib import reify
from .functions import length_of

def _to_array(value):
    """As a convenience, turn Python lists and tuples into NumPy arrays."""
    if isinstance(value, (tuple, list)):
        return array(value)
    elif isinstance(value, (float, int)):
        return np.float64(value)
    else:
        return value

class UnpackingError(Exception):
    """You cannot iterate directly over a Skyfield measurement object."""

class Distance(object):
    """A distance, stored internally as au and available in other units.

    You can initialize a ``Distance`` by providing a single float or a
    float array as either an ``au=``, ``km=``, or ``m=`` parameter.

    You can access the magnitude of the distance with its three
    attributes ``.au``, ``.km``, and ``.m``.  By default a distance
    prints itself in astronomical units (au), but you can take control
    of the formatting and choice of units yourself using standard Python
    numeric formatting:

    >>> d = Distance(au=1)
    >>> print(d)
    1.0 au
    >>> print('{:.2f} km'.format(d.km))
    149597870.70 km

    """
    _warned = False

    def __init__(self, au=None, km=None, m=None):
        if au is not None:
            self.au = _to_array(au)
        elif km is not None:
            self.km = _to_array(km)
            self.au = km / AU_KM
        elif m is not None:
            self.m = _to_array(m)
            self.au = m / AU_M
        else:
            raise ValueError('to construct a Distance provide au, km, or m')

    @classmethod
    def from_au(cls, au):
        self = cls.__new__(cls)
        self.au = _to_array(au)
        return self

    @reify
    def km(self):
        return self.au * AU_KM

    @reify
    def m(self):
        return self.au * AU_M

    @reify
    def AU(self):
        if not Distance._warned:
            print('WARNING: the IAU has renamed the astronomical unit to'
                  ' lowercase "au" so Skyfield will soon remove uppercase'
                  ' "AU" from Distance objects', file=sys.stdout)
            Distance._warned = True
        return self.au

    def __str__(self):
        n = self.au
        return ('{0} au' if getattr(n, 'shape', 0) else '{0:.6} au').format(n)

    def __repr__(self):
        return '<{0} {1}>'.format(type(self).__name__, self)

    def __iter__(self):
        cn = self.__class__.__name__
        raise UnpackingError(_iter_message.format(
            cls=cn, object=cn.lower(), values='x, y, z',
            attr1='au', attr2='km',
        ))

    def length(self):
        """Compute the length when this is an x,y,z vector.

        The Euclidean vector length of this vector is returned as a new
        :class:`~skyfield.units.Distance` object.

        >>> from skyfield.api import Distance
        >>> d = Distance(au=[1, 1, 0])
        >>> d.length()
        <Distance 1.41421 au>

        """
        return Distance(au=length_of(self.au))

    def light_seconds(self):
        """Return the length of this vector in light seconds."""
        return self.m / C

    def to(self, unit):
        """Convert this distance to the given AstroPy unit."""
        from astropy.units import au
        return (self.au * au).to(unit)

class Velocity(object):
    """A velocity, stored internally as au/day and available in other units.

    You can initialize a ``Velocity`` by providing a float or float
    array to its ``au_per_d=`` parameter.

    """
    _warned = False

    def __init__(self, au_per_d=None, km_per_s=None):
        if km_per_s is not None:
            au_per_d = km_per_s * DAY_S / AU_KM
        self.au_per_d = _to_array(au_per_d)

    @reify
    def km_per_s(self):
        return self.au_per_d * AU_KM / DAY_S

    @reify
    def m_per_s(self):
        return self.au_per_d * AU_M / DAY_S

    @reify
    def AU_per_d(self):
        if not Velocity._warned:
            print('WARNING: the IAU has renamed the astronomical unit to'
                  ' lowercase "au" so Skyfield will soon remove'
                  ' "AU_per_day" in favor of "au_per_day"',
                  file=sys.stdout)
            Velocity._warned = True
        return self.au_per_d

    def __str__(self):
        n = self.au_per_d
        fmt = '{0} au/day' if getattr(n, 'shape', 0) else '{0:.6} au/day'
        return fmt.format(n)

    def __repr__(self):
        return '<{0} {1}>'.format(type(self).__name__, self)

    def __iter__(self):
        cn = self.__class__.__name__
        raise UnpackingError(_iter_message.format(
            cls=cn, object=cn.lower(), values='xdot, ydot, zdot',
            attr1='au_per_d', attr2='km_per_s',
        ))

    def to(self, unit):
        """Convert this velocity to the given AstroPy unit."""
        from astropy.units import au, d
        return (self.au_per_d * au / d).to(unit)

_iter_message = """\
cannot directly unpack a {cls} into several values

To unpack a {cls} into three components, you need to ask for its
value in specific units through an attribute or method:

    {values} = {object}.{attr1}
    {values} = {object}.{attr2}
    {values} = {object}.to(astropy_unit)
"""

# Angle units.

_instantiation_instructions = """to instantiate an Angle, try one of:

Angle(angle=another_angle)
Angle(radians=value)
Angle(degrees=value)
Angle(hours=value)

where `value` can be either a Python float, a list of Python floats,
or a NumPy array of floats"""

class Angle(object):

    def __init__(self, angle=None, radians=None, degrees=None, hours=None,
                 preference=None, signed=False):

        if angle is not None:
            if not isinstance(angle, Angle):
                raise ValueError(_instantiation_instructions)
            self.radians = angle.radians
        elif radians is not None:
            self.radians = _to_array(radians)
        elif degrees is not None:
            self._degrees = degrees = _to_array(_unsexagesimalize(degrees))
            self.radians = degrees / 360.0 * tau
        elif hours is not None:
            self._hours = hours = _to_array(_unsexagesimalize(hours))
            self.radians = hours / 24.0 * tau

        self.preference = (preference if preference is not None
                           else 'hours' if hours is not None
                           else 'degrees')
        self.signed = signed

    @classmethod
    def from_degrees(cls, degrees, signed=False):
        degrees = _to_array(_unsexagesimalize(degrees))
        self = cls.__new__(cls)
        self.degrees = degrees
        self.radians = degrees / 360.0 * tau
        self.preference = 'degrees'
        self.signed = signed
        return self

    @reify
    def _hours(self):
        return self.radians * 24.0 / tau

    @reify
    def _degrees(self):
        return self.radians * 360.0 / tau

    @reify
    def hours(self):
        if self.preference != 'hours':
            raise WrongUnitError('hours')
        return self._hours

    @reify
    def degrees(self):
        if self.preference != 'degrees':
            raise WrongUnitError('degrees')
        return self._degrees

    def arcminutes(self):
        """Return the angle in arcminutes."""
        return self._degrees * 60.0

    def arcseconds(self):
        """Return the angle in arcseconds."""
        return self._degrees * 3600.0

    def mas(self):
        """Return the angle in milliarcseconds."""
        return self._degrees * 3600000.0

    def __str__(self):
        if self.radians.size == 0:
            return 'Angle []'
        return self.dstr() if self.preference == 'degrees' else self.hstr()

    def __repr__(self):
        if self.radians.size == 0:
            return '<{0} []>'.format(type(self).__name__)
        else:
            return '<{0} {1}>'.format(type(self).__name__, self)

    def __iter__(self):
        raise ValueError(
            '''choose a specific Angle unit to iterate over

Instead of iterating over this Angle object, try iterating over one of
its unit-specific arrays like .degrees, .hours, or .radians, or else over
the output of one of its methods like .hstr(), .dstr(), .arcminutes(),
.arcseconds(), or .mas().  For all of the possibilities see:
https://rhodesmill.org/skyfield/api-units.html#skyfield.units.Angle''')

    def hms(self, warn=True):
        """Convert to a tuple (hours, minutes, seconds).

        All three quantities will have the same sign as the angle itself.

        """
        if warn and self.preference != 'hours':
            raise WrongUnitError('hms')
        sign, units, minutes, seconds = _sexagesimalize_to_float(self._hours)
        return sign * units, sign * minutes, sign * seconds

    def signed_hms(self, warn=True):
        """Convert to a tuple (sign, hours, minutes, seconds).

        The ``sign`` will be either +1 or -1, and the other quantities
        will all be positive.

        """
        if warn and self.preference != 'hours':
            raise WrongUnitError('signed_hms')
        return _sexagesimalize_to_float(self._hours)

    def hstr(self, places=2, warn=True):
        """Convert to a string like ``12h 07m 30.00s``."""
        if warn and self.preference != 'hours':
            raise WrongUnitError('hstr')
        if self.radians.size == 0:
            return '<Angle []>'
        hours = self._hours
        shape = getattr(hours, 'shape', ())
        if shape and shape != (1,):
            return "{0} values from {1} to {2}".format(
                len(hours),
                _hstr(hours[0], places),
                _hstr(hours[-1], places),
                )
        return _hstr(hours, places)

    def dms(self, warn=True):
        """Convert to a tuple (degrees, minutes, seconds).

        All three quantities will have the same sign as the angle itself.

        """
        if warn and self.preference != 'degrees':
            raise WrongUnitError('dms')
        sign, units, minutes, seconds = _sexagesimalize_to_float(self._degrees)
        return sign * units, sign * minutes, sign * seconds

    def signed_dms(self, warn=True):
        """Convert to a tuple (sign, degrees, minutes, seconds).

        The ``sign`` will be either +1 or -1, and the other quantities
        will all be positive.

        """
        if warn and self.preference != 'degrees':
            raise WrongUnitError('signed_dms')
        return _sexagesimalize_to_float(self._degrees)

    def dstr(self, places=1, warn=True):
        """Convert to a string like ``181deg 52\' 30.0"``."""
        if warn and self.preference != 'degrees':
            raise WrongUnitError('dstr')
        if self.radians.size == 0:
            return '<Angle []>'
        degrees = self._degrees
        signed = self.signed
        shape = getattr(degrees, 'shape', ())
        if shape and shape != (1,):
            return "{0} values from {1} to {2}".format(
                len(degrees),
                _dstr(degrees[0], places, signed),
                _dstr(degrees[-1], places, signed),
                )
        return _dstr(degrees, places, signed)

    def to(self, unit):
        """Convert this angle to the given AstroPy unit."""
        from astropy.units import rad
        return (self.radians * rad).to(unit)

        # Or should this do:
        from astropy.coordinates import Angle
        from astropy.units import rad
        return Angle(self.radians, rad).to(unit)

class WrongUnitError(ValueError):

    def __init__(self, name):
        unit = 'hours' if (name.startswith('h') or '_h' in name) else 'degrees'
        usual = 'hours' if (unit == 'degrees') else 'degrees'
        message = ('this angle is usually expressed in {0}, not {1};'
                   ' if you want to use {1} anyway,'.format(usual, unit))
        if name == unit:
            message += ' then please use the attribute _{0}'.format(unit)
        else:
            message += ' then call {0}() with warn=False'.format(name)
        self.args = (message,)

def _sexagesimalize_to_float(value):
    """Decompose `value` into units, minutes, and seconds.

    Note that this routine is not appropriate for displaying a value,
    because rounding to the smallest digit of display is necessary
    before showing a value to the user.  Use `_sexagesimalize_to_int()`
    for data being displayed to the user.

    This routine simply decomposes the floating point `value` into a
    sign (+1.0 or -1.0), units, minutes, and seconds, returning the
    result in a four-element tuple.

    >>> _sexagesimalize_to_float(12.05125)
    (1.0, 12.0, 3.0, 4.5)
    >>> _sexagesimalize_to_float(-12.05125)
    (-1.0, 12.0, 3.0, 4.5)

    """
    sign = np.sign(value)
    n = abs(value)
    minutes, seconds = divmod(n * 3600.0, 60.0)
    units, minutes = divmod(minutes, 60.0)
    return sign, units, minutes, seconds

def _sexagesimalize_to_int(value, places=0):
    """Decompose `value` into units, minutes, seconds, and second fractions.

    This routine prepares a value for sexagesimal display, with its
    seconds fraction expressed as an integer with `places` digits.  The
    result is a tuple of five integers:

    ``(sign [either +1 or -1], units, minutes, seconds, second_fractions)``

    The integers are properly rounded per astronomical convention so
    that, for example, given ``places=3`` the result tuple ``(1, 11, 22,
    33, 444)`` means that the input was closer to 11u 22' 33.444" than
    to either 33.443" or 33.445" in its value.

    """
    sign = int(np.sign(value))
    value = abs(value)
    power = 10 ** places
    n = int(7200 * power * value + 1) // 2
    n, fraction = divmod(n, power)
    n, seconds = divmod(n, 60)
    n, minutes = divmod(n, 60)
    return sign, n, minutes, seconds, fraction

def _hstr(hours, places=2):
    """Convert floating point `hours` into a sexagesimal string.

    >>> _hstr(12.125)
    '12h 07m 30.00s'
    >>> _hstr(12.125, places=4)
    '12h 07m 30.0000s'
    >>> _hstr(float('nan'))
    'nan'

    """
    if isnan(hours):
        return 'nan'
    sgn, h, m, s, etc = _sexagesimalize_to_int(hours, places)
    sign = '-' if sgn < 0.0 else ''
    return '%s%02dh %02dm %02d.%0*ds' % (sign, h, m, s, places, etc)

def _dstr(degrees, places=1, signed=False):
    r"""Convert floating point `degrees` into a sexagesimal string.

    >>> _dstr(181.875)
    '181deg 52\' 30.0"'
    >>> _dstr(181.875, places=3)
    '181deg 52\' 30.000"'
    >>> _dstr(181.875, signed=True)
    '+181deg 52\' 30.0"'
    >>> _dstr(float('nan'))
    'nan'

    """
    if isnan(degrees):
        return 'nan'
    sgn, d, m, s, etc = _sexagesimalize_to_int(degrees, places)
    sign = '-' if sgn < 0.0 else '+' if signed else ''
    return '%s%02ddeg %02d\' %02d.%0*d"' % (sign, d, m, s, places, etc)

def wms(whole, minutes=0.0, seconds=0.0):
    """Return a quantity expressed with 1/60 minutes and 1/3600 seconds."""
    return (whole
            + copysign(minutes, whole) / 60.0
            + copysign(seconds, whole) / 3600.0)

def _unsexagesimalize(value):
    """Return `value` after interpreting a (units, minutes, seconds) tuple.

    When `value` is not a tuple, it is simply returned.

    >>> _unsexagesimalize(3.25)
    3.25

    An input tuple is interpreted as units, minutes, and seconds.  Note
    that only the sign of `units` is significant!  So all of the
    following tuples convert into exactly the same value:

    >>> '%f' % _unsexagesimalize((-1, 2, 3))
    '-1.034167'
    >>> '%f' % _unsexagesimalize((-1, -2, 3))
    '-1.034167'
    >>> '%f' % _unsexagesimalize((-1, -2, -3))
    '-1.034167'

    """
    if isinstance(value, tuple):
        components = iter(value)
        value = next(components)
        factor = 1.0
        for component in components:
            factor *= 60.0
            value += copysign(component, value) / factor
    return value

def _interpret_angle(name, angle_object, angle_float, unit='degrees'):
    """Return an angle in radians from one of two arguments.

    It is common for Skyfield routines to accept both an argument like
    `alt` that takes an Angle object as well as an `alt_degrees` that
    can be given a bare float or a sexagesimal tuple.  A pair of such
    arguments can be passed to this routine for interpretation.

    """
    if angle_object is not None:
        if isinstance(angle_object, Angle):
            return angle_object.radians
    elif angle_float is not None:
        return _unsexagesimalize(angle_float) / 360.0 * tau
    raise ValueError('you must either provide the {0}= parameter with'
                     ' an Angle argument or supply the {0}_{1}= parameter'
                     ' with a numeric argument'.format(name, unit))

def _ltude(value, name, psuffix, nsuffix):
    # Support for old deprecated Topos argument interpretation.
    if not isinstance(value, str):
        return _unsexagesimalize(value)
    value = value.strip().upper()
    if value.endswith(psuffix):
        sign = +1.0
    elif value.endswith(nsuffix):
        sign = -1.0
    else:
        raise ValueError('your {0} string {1!r} does not end with either {2!r}'
                         ' or {3!r}'.format(name, value, psuffix, nsuffix))
    try:
        value = float(value[:-1])
    except ValueError:
        raise ValueError('your {0} string {1!r} cannot be parsed as a floating'
                         ' point number'.format(name, value))
    return sign * value
