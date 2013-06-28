"""Python classes that represent various classes of star."""

from .timescales import T0

class Star(object):

    def __init__(self, ra, dec,
                 pm_ra=0.0, pm_dec=0.0, parallax=0.0, rad_vel=0.0):
        self.ra = ra
        self.dec = dec
        self.pm_ra = pm_ra
        self.pm_dec = pm_dec
        self.parallax = parallax
        self.rad_vel = rad_vel

    def __call__(self, jd):
        position, velocity = 'starvectors(self)'  # pos1
        dt = 'd_light'(position, 'observer')
        position = 'proper_motion'(T0, position, velocity, (jd + dt))
        position, t_light = 'bary2obs(position, observer)'
        return position
