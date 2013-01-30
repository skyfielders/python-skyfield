
============
Introduction
============


>>> import sgp4


>>> from ephem.planets import Ephemeris
>>> eph = Ephemeris()
>>> earth, mars = eph.earth, eph.mars
>>> print earth(2414993.5).observe(mars).astrometric()


>>> import numpy as np
>>> t0 = 2414993.5
>>> t = np.arange(t0, t0 + 5, 1.0)
>>> print earth(t).observe(mars).astrometric()

>>> print earth(t).observe(mars).apparent()

        ggr = coordinates.Topos('75 W', '45 N', 0.0,
                                temperature=10.0, pressure=1010.0)
        ggr.earth = self.e.earth
        ggr.ephemeris = self.e
        delta_t = 0

        for t, name in product((T0, TA, TB), planets_to_test):
            obj = c.make_object(0, planet_codes[name], b'planet', None)
            ra, dec, dis = c.topo_planet(t, delta_t, obj, position)

            planet = getattr(self.e, name)
            g = ggr(t).observe(planet).apparent()
