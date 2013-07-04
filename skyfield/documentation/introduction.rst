
============
Introduction
============


>>> import sgp4
>>> import numpy as np

>>> from skyfield.planets import Ephemeris
>>> eph = Ephemeris()
>>> earth, mars = eph.earth, eph.mars

>>> e = earth(2414993.5)
>>> print e
<ICRS position x,y,z AU and velocity xdot,ydot,zdot AU/day at date jd>
>>> print e.x
0.273833256524
>>> print e.y
0.874908499106
>>> print e.z
0.37944053759

>>> print e.xdot
-0.0168396519412
>>> print e.ydot
0.00427785895254
>>> print e.zdot
0.00185493319076

>>> dates = np.array([2414993.5, 2414994.5])
>>> e = earth(dates)
>>> e.x
array([ 0.27383326,  0.25695284])

huh, how do we measure distance? abs does not work >>> print
  mars(2414993.5).position - earth(2414993.5).position

>>> emi = mars(2414993.5) - earth(2414993.5)
>>> print emi.x, emi.y, emi.z
-0.220661979709 -2.18621208383 -0.982457068395

>>> print mars.observe_from(earth(2414993.5))
<GCRS position x,y,z AU and velocity xdot,ydot,zdot AU/day at date jd>

>>> print mars.observe_from(earth(2414993.5)).astrometric()
<Astrometric position RA=4.61170587399 dec=-0.420447658663>
>>> print mars.observe_from(earth(2414993.5)).apparent()
<Apparent position RA=[ 4.58500039] dec=[-0.41933084]>


>>> a = mars.observe_from(earth(2414993.5)).astrometric()
>>> a
<Astrometric position RA=4.61170587399 dec=-0.420447658663>
>>> a.ra
HourAngle(4.611705873992877)
>>> a.ra.hms()
(1.0, 17.0, 36.0, 55.507904515793314)
>>> a.ra.hstr()
'17h 36m 55.5079045158s'
>>> a.dec.dms()
(-1.0, 24.0, 5.0, 23.554851165623347)
>>> a.dec.dstr()
'-24deg 5m 23.5548511656s'

>>> a.ra.degrees()
Traceback (most recent call last):
  ...
WrongUnitError: This angle is usually expressed in hours, not degrees; if you want to express it in degrees anyway, use degrees_anyway()

>>> a.dec.hms()
Traceback (most recent call last):
  ...
WrongUnitError: This angle is usually expressed in degrees, not hours; if you want to express it in hours anyway, use hms_anyway()

repr(a.dec.dpretty())
-24
-24°5´23´´.5548511656

>>> a = mars.observe_from(earth(np.array([2414993.5, 2414994.5]))).astrometric()
>>> a
<Astrometric position RA=[ 4.61170587  4.62603456] dec=[-0.42044766 -0.42106794]>


>>> import numpy as np
>>> t0 = 2414993.5
>>> t = np.arange(t0, t0 + 5, 1.0)
>>> print mars.observe_from(earth(t)).astrometric()
<Astrometric position RA=[ 4.61170587  4.62603456  4.64038665  4.65476093  4.66915621] dec=[-0.42044766 -0.42106794 -0.42161316 -0.42208295 -0.42247693]>

>>> print mars.observe_from(earth(t)).apparent()
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
            g = planet.observe_from(ggr(t)).apparent()


DE405  52.1 MB  1600–2200 (May 1997)
DE406 170.0 MB -3000–3000 (May 1997)
DE421  13.0 MB  1900–2050 (February 2008)
DE422 519.6 MB -3000–3000 (September 2009)
DE423  34.6 MB  1800–2200 (February 2010)

