
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


DE405  52.1 MB  1600–2200 (May 1997)
DE406 170.0 MB -3000–3000 (May 1997)
DE421  13.0 MB  1900–2050 (February 2008)
DE422 519.6 MB -3000–3000 (September 2009)
DE423  34.6 MB  1800–2200 (February 2010)

