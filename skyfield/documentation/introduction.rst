
============
Introduction
============


>>> import sgp4


>>> from skyfield.planets import Ephemeris
>>> eph = Ephemeris()
>>> earth, mars = eph.earth, eph.mars
>>> print earth(2414993.5)
<ICRS position x=[ 0.27383326] y=[ 0.8749085] z=[ 0.37944054]>
>>> print earth(2414993.5).observe(mars)
<GCRS position x=[-0.22086377] y=[-2.1862353] z=[-0.98246221]>
>>> print earth(2414993.5).observe(mars).astrometric()
<Astrometric position RA=[ 4.61170587] dec=[-0.42044766]>

>>> import numpy as np
>>> t0 = 2414993.5
>>> t = np.arange(t0, t0 + 5, 1.0)
>>> print earth(t).observe(mars).astrometric()
<Astrometric position RA=[ 4.61170587  4.62603456  4.64038665  4.65476093  4.66915621] dec=[-0.42044766 -0.42106794 -0.42161316 -0.42208295 -0.42247693]>

>>> print earth(t).observe(mars).apparent()
<Apparent position RA=[ 4.58500039  4.59931574  4.61365598  4.62801986  4.64240623] dec=[-0.41933084 -0.42008982 -0.42077423 -0.42138359 -0.42191745]>


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

