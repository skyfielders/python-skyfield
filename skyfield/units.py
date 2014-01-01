from .constants import AU_KM

class Distance(object):

    def __init__(self, value_AU):
        self.AU = value_AU

    def __getattr__(self, name):
        if name == 'km':
            self.km = self.AU * AU_KM
            return self.km
        return super(Distance, self).__getattr__(name)
