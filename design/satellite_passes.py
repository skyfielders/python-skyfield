import numpy as np
from skyfield import almanac, api
from skyfield.nutationlib import iau2000b

ts = api.load.timescale(builtin=True)
eph = api.load('de421.bsp')
sat = api.EarthSatellite(
    '1 25544U 98067A   20150.54009403  .00000886  00000-0  23936-4 0  9992',
    '2 25544  51.6447  80.7632 0002782   3.7438 167.9307 15.49399736229124',
)

topos = api.Topos('18.3368 N', '64.7281 W')
t0 = ts.utc(2020, 6, 5)
t1 = ts.utc(2020, 6, 6)

half_second = 0.5 / 24.0 / 3600.0

from time import time

T0 = time()
t, y = sat.find_events(topos, t0, t1, altitude_degrees=10.0)
print(time() - T0, 'seconds to find events:')

for ti, yi in zip(t, y):
    print(' ', ti.utc_jpl(), yi)

offset = 0
assert list(y[offset:offset+3]) == [0, 1, 2]

pass_start, pass_end = t[offset+0], t[offset+2]

def f(t):
    length = t.tt.shape[0]
    t._nutation_angles = iau2000b(t.tt)
    return sat.at(t).is_sunlit(eph)

f.rough_period = 0.001

T0 = time()
t2, y2 = almanac.find_discrete(pass_start, pass_end, f, epsilon=half_second)
print(time() - T0, 'seconds to find moment of entering shadow:')

for ti, yi in zip(t2, y2):
    print(' ', ti.utc_jpl(), 'entered sunlight' if yi else 'entered shadow')

one_second = 1.0 / 24.0 / 3600.0
tt = np.arange(pass_start.tt, pass_end.tt, one_second)

T0 = time()
t = ts.tt_jd(tt)
t._nutation_angles = iau2000b(t.tt)
satpos = (sat - topos).at(t)
is_sunlit = satpos.is_sunlit(eph)
DT = time() - T0

print(DT, 'to compute is_sunlit() for every second of {}-second pass'
      .format(len(tt)))

T0 = time()
alt, az, distance = satpos.altaz()
DT = time() - T0

print(DT, 'to re-use the positions to compute altitude and azimuth')

print('''
On my laptop, this script shows that simply computing the positions and
the sunlit-ness of the ISS for each of the 381 seconds of this pass
takes roughly the same amount of time as mounting a full search for the
moment it passes into shadow.  But the benefit is far greater, because
with almost no additional expense all altitudes and azimuths can also
be computed.
''')
