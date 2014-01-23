from .constants import AU_KM, DAY_S

class Distance(object):
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

    def to(self, unit):
        from astropy.units import AU
        return (self.AU * AU).to(unit)

class Velocity(object):
    def __init__(self, AU_per_d):
        self.AU_per_d = AU_per_d

    def __getattr__(self, name):
        if name == 'km_per_s':
            self.km_per_s = self.AU_per_d * AU_KM / DAY_S
            return self.km_per_s
        raise AttributeError('no attribute named %r' % (name,))

    def to(self, unit):
        from astropy.units import AU, d
        return (self.AU_per_d * AU / d).to(unit)
