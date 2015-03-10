"""Accuracy tests against data pulled from HORIZONS."""

from numpy import max
from skyfield import api
from skyfield.constants import AU_M
from skyfield.io import download
from skyfield.jpllib import Kernel

one_second = 1.0 / 24.0 / 60.0 / 60.0
arcsecond = 1.0 / 60.0 / 60.0
ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
meter = 1.0 / AU_M

base = 'http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk'
de430_url = base + '/planets/de430.bsp'
de431_url = base + '/planets/de431.bsp'
jup310_url = base + '/satellites/jup310.bsp'

def compare(value, expected_value, epsilon):
    if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
        assert max(abs(value - expected_value)) <= epsilon
    else:
        assert abs(value - expected_value) <= epsilon

def test_jupiter1():
    astrometric = api.sun(utc=(1980, 1, 1, 0, 0)).observe(api.jupiter)
    hlat, hlon, d = astrometric.ecliptic_latlon()
    compare(hlat.degrees, 1.013, 0.001)
    compare(hlon.degrees, 151.3229, 0.001)

def test_callisto_geometry():
    k = Kernel(download(jup310_url))
    a = k.earth.geometry_of(k.callisto).at(tdb=2471184.5)
    compare(a.position.au,
      [-4.884815926454119E+00, -3.705745549073268E+00, -1.493487818022234E+00],
      0.001 * meter)
    compare(a.velocity.au_per_d,
      [9.604665478763035E-03, -1.552997751083403E-02, -6.678445860769302E-03],
      0.000001 * meter)

def test_callisto_astrometric():
    k = Kernel(download(jup310_url))
    a = k.earth.observe(k.callisto).at(utc=(2053, 10, 9))
    ra, dec, distance = a.radec()
    compare(ra._degrees, 217.1839292, 0.001 * arcsecond)
    compare(dec.degrees, -13.6892791, 0.001 * arcsecond)
    compare(distance.au, 6.31079291776184, 0.1 * meter)

def other():
    jd = api.JulianDate(tdb=(2015, 3, 2), delta_t=67.185390)
    bos = api.earth.topos(
        '42.3567000 N', '288.943100 E') #, elevation_m=92.0)

    a = api.earth(jd) - bos(jd) #bos(jd).observe(api.earth)
    v = a.position.au
    print('usno', [repr(n) for n in v / meter])

    from skyfield.functions import length_of
    print('r = ', length_of(v) / meter)
    from ..constants import C
    print('r(c) = ', length_of(v) / meter / C)
    from skyfield.constants import ERAD
    print('R = ', ERAD)

def off_test_boston_geometry():
    print('-------------------------------')
    other()
    print('-------')
    k = Kernel(download(de430_url))
    #boston = api.Topos((42, 21, 24.1), (-71, 3, 24.8))
    boston = api.Topos('42.3567000 N', '288.943100 E', elevation_m=0.0)#1.014E-9)
    #print(boston.longitude.degrees+360.0)
    #print(boston.observe(k.moon))
    #print(boston.observe(k.earth))
    jd = api.JulianDate(tdb=(2015, 3, 2), delta_t=67.185390)
    #a = boston.observe(k.earth).at(tdb=(2015, 3, 2))
    a = boston.geometry_of(k.earth).at(jd) #tdb=(2015, 3, 2, 0, -1, -7.71))
    #print(a.jd.tdb)
    print('got ', [repr(n) for n in (a.position.au) / meter])
    from numpy import array
    w = array([-1.179627735416499E-06, -3.153207525942794E-05, -2.857611787776292E-05])
    print('want', [repr(n) for n in (w) / meter])
    print('diff', [repr(n) for n in (a.position.au - w) / meter])
    #print('at() light time:', a.lighttime)
    from skyfield.functions import length_of
    print(length_of(a.position.au) / meter)
    print(length_of(w) / meter)
    #
    compare(a.position.au,
      [-1.179627735416499E-06, -3.153207525942794E-05, -2.857611787776292E-05],
      0.001 * meter)

def off_test_moon_from_boston_geometry():
    k = Kernel(download(de430_url))
    boston = api.Topos((42, 21, 24.1), (-71, 3, 24.8))
    boston = api.Topos('42.3567000 N', '288.943100 E')
    #print(boston.longitude.degrees+360.0)
    #print(boston.observe(k.moon))
    jd = api.JulianDate(tdb=(2015, 3, 2), delta_t=67.185390)
    #a = boston.observe(k.moon).at(jd)
    a = boston.geometry_of(k.moon).at(jd)
    diff = a.position.au - [-1.341501206552443E-03, 2.190483327459023E-03, 6.839177007993498E-04]
    print('diff (au)', [repr(n) for n in diff])
    print('diff (m) ', [repr(n) for n in diff * AU_M])
    #
    compare(a.position.au,
      [-1.341501206552443E-03, 2.190483327459023E-03, 6.839177007993498E-04],
      0.001 * meter)

def off_test_moon_from_boston_astrometric():
    k = Kernel(download(de430_url))
    boston = api.Topos((42, 21, 24.1), (-71, 3, 24.8))
    boston = api.Topos('42.3567000 N', '288.943100 E')
    #print(boston.longitude.degrees+360.0)
    #print(boston.observe(k.moon))
    a = boston.observe(k.moon).at(utc=(2015, 3, 2))#.apparent()
    #a = boston.geometry_of(k.moon).at(utc=(2015, 3, 2))
    ra, dec, distance = a.radec()
    print(ra._degrees, ra._degrees - 121.4863918)
    print(dec._degrees, dec.degrees - 14.9096363)
    compare(ra._degrees, 121.4863918, 0.001 * arcsecond)
    compare(dec.degrees, 14.9096363, 0.001 * arcsecond)
    compare(distance.au, 0.00265821944771, 0.1 * meter)
    # azi el 123.7693  50.3373
