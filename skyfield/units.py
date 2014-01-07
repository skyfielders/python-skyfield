from .constants import AU_KM, DAY_S

class Distance(object):
    def __init__(self, value_AU):
        self.AU = value_AU

    def __getattr__(self, name):
        if name == 'km':
            self.km = self.AU * AU_KM
            return self.km
        raise AttributeError('no attribute named %r' % (name,))

class Velocity(object):
    def __init__(self, value_AU_per_d):
        self.AU_per_d = value_AU_per_d

    def __getattr__(self, name):
        if name == 'km_per_s':
            self.km_per_s = self.AU_per_d * AU_KM / DAY_S
            return self.km_per_s
        raise AttributeError('no attribute named %r' % (name,))
