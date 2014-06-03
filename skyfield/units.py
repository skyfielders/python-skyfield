"""Simple distance, velocity, and angle support for Skyfield.

"""
import numpy as np
from .constants import AU_KM, DAY_S, tau

# Distance and velocity.

class UnpackingError(Exception):
    """You cannot iterate directly over a Skyfield measurement object."""

class Distance(object):
    """A distance, stored internally as AU and available in other units.

    You can initialize a ``Distance`` by providing a single float or a
    float array as either an ``AU=`` parameter or a ``km=`` parameter
    when building a ``Distance`` object.

    """
    def __init__(self, AU=None, km=None):
        if AU is not None:
            self.AU = AU
        elif km is not None:
            self.km = km
            self.AU = km / AU_KM

    def __getattr__(self, name):
        if name == 'km':
            self.km = self.AU * AU_KM
            return self.km
        raise AttributeError('no attribute named %r' % (name,))

    def __str__(self):
        return '%s AU' % self.AU

    def __iter__(self):
        raise UnpackingError(_iter_message % {
            'class': self.__class__.__name__, 'values': 'x, y, z',
            'attr1': 'AU', 'attr2': 'km'})

    def to(self, unit):
        """Return this distance in the given AstroPy units."""
        from astropy.units import AU
        return (self.AU * AU).to(unit)

class Velocity(object):
    """A velocity, stored internally as AU/day and available in other units.

    You can initialize a ``Velocity`` by providing a single float or a
    float array as either an ``AU_per_d=`` parameter.

    """
    def __init__(self, AU_per_d):
        self.AU_per_d = AU_per_d

    def __getattr__(self, name):
        if name == 'km_per_s':
            self.km_per_s = self.AU_per_d * AU_KM / DAY_S
            return self.km_per_s
        raise AttributeError('no attribute named %r' % (name,))

    def __str__(self):
        return '%s AU/day' % self.AU_per_d

    def __iter__(self):
        raise UnpackingError(_iter_message % {
            'class': self.__class__.__name__, 'values': 'xdot, ydot, zdot',
            'attr1': 'AU_per_d', 'attr2': 'km_per_s'})

    def to(self, unit):
        """Return this velocity in the given AstroPy units."""
        from astropy.units import AU, d
        return (self.AU_per_d * AU / d).to(unit)

_iter_message = """\
cannot directly unpack a %(class)s into several values

To unpack a %(class)s into three components, you need to ask for its
value in specific units through an attribute or method:

    %(values)s = velocity.%(attr1)s
    %(values)s = velocity.%(attr2)s
    %(values)s = velocity.to(astropy_unit)
"""

# Angle units.

_to_degrees = 360.0 / tau
_from_degrees = tau / 360.0

_to_hours = 24.0 / tau
_from_hours = tau / 24.0

_instantiation_instructions = """to instantiate an Angle, try one of:

Angle(angle=another_angle)
Angle(radians=value)
Angle(degrees=value)
Angle(hours=value)

where `value` can be either a Python float or a NumPy array of floats"""

class BaseAngle(object):

    _unary_plus = False

    def __init__(self, angle=None, radians=None, degrees=None, hours=None):
        if angle is not None:
            if not isinstance(angle, BaseAngle):
                raise ValueError(_instantiation_instructions)
            self._radians = angle._radians
        elif radians is not None:
            self._radians = radians
        elif degrees is not None:
            self._radians = degrees * _from_degrees
        elif hours is not None:
            self._radians = hours * _from_hours

    def __format__(self, format_spec):
        return self.dstr()

    def radians(self):
        return self._radians

    def hours(self):
        return self._radians * _to_hours

    def hms(self):
        return _sexagesimalize(self.hours())

    def hstr(self, places=2, plus=False):
        sgn, h, m, s, etc = _sexagesimalize(self.hours(), places)
        sign = '-' if sgn < 0.0 else '+' if (plus or self._unary_plus) else ''
        return '%s%02dh %02dm %02d.%0*ds' % (sign, h, m, s, places, etc)

    def degrees(self):
        return self._radians * _to_degrees

    def dms(self):
        return _sexagesimalize(self.degrees())

    def dstr(self, places=1, plus=False):
        sgn, d, m, s, etc = _sexagesimalize(self.degrees(), places)
        sign = '-' if sgn < 0.0 else '+' if (plus or self._unary_plus) else ''
        return '%s%02ddeg %02d\' %02d.%0*d"' % (sign, d, m, s, places, etc)

    hours_anyway = hours
    hms_anyway = hms
    hstr_anyway = hstr

    degrees_anyway = degrees
    dms_anyway = dms
    dstr_anyway = dstr

class WrongUnitError(ValueError):

    def __init__(self, name, unit):
        usual = 'hours' if (unit == 'degrees') else 'degrees'
        self.args = ('This angle is usually expressed in {}, not {};'
                     ' if you want to express it in {} anyway, use'
                     ' {}_anyway()'.format(usual, unit, unit, name),)

class Angle(BaseAngle):

    __str__ = BaseAngle.dstr

    # Protect naive users from accidentally calling hour methods.

    def hours(self):
        raise WrongUnitError('hours', 'hours')

    def hms(self):
        raise WrongUnitError('hms', 'hours')

    def hstr(self):
        raise WrongUnitError('hstr', 'hours')

class SignedAngle(Angle):
    """An Angle that prints a unary ``'+'`` when positive."""

    _unary_plus = True

class HourAngle(BaseAngle):

    __str__ = BaseAngle.hstr

    # Protect naive users from accidentally calling degree methods.

    def degrees(self):
        raise WrongUnitError('degrees', 'degrees')

    def dms(self):
        raise WrongUnitError('dms', 'degrees')

    def dstr(self):
        raise WrongUnitError('dstr', 'degrees')

def _sexagesimalize(value, places=0):
    sign = int(np.sign(value))
    value = np.absolute(value)
    power = 10 ** places
    n = int(7200 * power * value + 1) // 2
    n, fraction = divmod(n, power)
    n, seconds = divmod(n, 60)
    n, minutes = divmod(n, 60)
    return sign, n, minutes, seconds, fraction

def interpret_longitude(value):
    split = getattr(value, 'split', None)
    if split is not None:
        pieces = split()
        degrees = float(pieces[0])
        if len(pieces) > 1 and pieces[1].lower() == 'w':
            degrees = - degrees
        return degrees / 360. * tau
    else:
        return value

def interpret_latitude(value):
    split = getattr(value, 'split', None)
    if split is not None:
        pieces = split()
        degrees = float(pieces[0])
        if len(pieces) > 1 and pieces[1].lower() == 's':
            degrees = - degrees
        return degrees / 360. * tau
    else:
        return value
