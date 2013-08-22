
==============
 Introduction
==============


>>> from skyfield.timescales import JulianDate
>>> jd = JulianDate(ut1=2456755.75)

>>> from skyfield.planets import earth, mars
>>> boston = earth.topos('71.0636 W', '42.3583 N')
>>> h = boston(jd).observe(mars).apparent().horizontal()
>>> print(h.alt.dstr())
40deg 3m 49.433s
>>> print(h.az.dstr())
201deg 58m 7.067s

fix >>> import sgp4
fix >>> import numpy as np

fix >>> import skyfield.angles

fix >>> e = earth(2414993.5)
fix >>> print(e)
<ICRS position x,y,z AU and velocity xdot,ydot,zdot AU/day at date jd>
fix >>> print(e.x)
0.273833256524
fix >>> print(e.y)
0.874908499106
fix >>> print(e.z)
0.37944053759

fix >>> print(e.xdot)
-0.0168396519412
fix >>> print(e.ydot)
0.00427785895254
fix >>> print(e.zdot)
0.00185493319076

fix >>> dates = np.array([2414993.5, 2414994.5])
fix >>> e = earth(dates)
fix >>> e.x
array([ 0.27383326,  0.25695284])

huh, how do we measure distance? abs does not work fix >>> print
  mars(2414993.5).position - earth(2414993.5).position

fix >>> emi = mars(2414993.5) - earth(2414993.5)
fix >>> print(emi.position)
[-0.22066198 -2.18621208 -0.98245707]

fix >>> print(mars.observe_from(earth(2414993.5)))
<GCRS position x,y,z AU and velocity xdot,ydot,zdot AU/day at date jd>

fix >>> print(mars.observe_from(earth(2414993.5)).astrometric())
<Astrometric position RA=HourAngle(4.611705873992877) dec=Angle(-0.4204476586629836)>
fix >>> print(mars.observe_from(earth(2414993.5)).apparent())
<Apparent position RA=HourAngle([ 4.58500039]) dec=Angle([-0.41933084])>


fix >>> a = mars.observe_from(earth(2414993.5)).astrometric()
fix >>> a
<Astrometric position RA=HourAngle(4.611705873992877) dec=Angle(-0.4204476586629836)>
fix >>> a.ra
HourAngle(4.611705873992877)
fix >>> a.ra.hms()
(1.0, 17.0, 36.0, 55.507904515793314)
fix >>> a.ra.hstr()
'17h 36m 55.508s'
fix >>> a.dec.dms()
(-1.0, 24.0, 5.0, 23.554851165623347)
fix >>> a.dec.dstr()
'-24deg 5m 23.555s'

fix >>> try:
...     a.ra.degrees()
... except skyfield.angles.WrongUnitError as e:
...     print(e)
This angle is usually expressed in hours, not degrees; if you want to express it in degrees anyway, use degrees_anyway()

fix >>> try:
...     a.dec.hms()
... except skyfield.angles.WrongUnitError as e:
...     print(e)
This angle is usually expressed in degrees, not hours; if you want to express it in hours anyway, use hms_anyway()

repr(a.dec.dpretty())
-24
-24°5´23´´.5548511656

fix >>> a = mars.observe_from(earth(np.array([2414993.5, 2414994.5]))).astrometric()
fix >>> a
<Astrometric position RA=HourAngle([ 4.61170587,  4.62603456]) dec=Angle([-0.42044766, -0.42106794])>


fix >>> import numpy as np
fix >>> t0 = 2414993.5
fix >>> t = np.arange(t0, t0 + 5, 1.0)
fix >>> print(mars.observe_from(earth(t)).astrometric())
<Astrometric position RA=HourAngle([ 4.61170587,  4.62603456,  4.64038665,  4.65476093,  4.66915621]) dec=Angle([-0.42044766, -0.42106794, -0.42161316, -0.42208295, -0.42247693])>

fix >>> print(mars.observe_from(earth(t)).apparent())
<Apparent position RA=HourAngle([ 4.58500039,  4.59931574,  4.61365598,  4.62801986,  4.64240623]) dec=Angle([-0.41933084, -0.42008982, -0.42077423, -0.42138359, -0.42191745])>


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

