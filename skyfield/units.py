from .constants import AU_KM, DAY_S

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
