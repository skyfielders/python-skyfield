#!/usr/bin/env python3

'''
Compute time: 0.163 seconds
Altitude vs horizon: min -5.176029 arcseconds,  max 2.819515 arcseconds
Worst undershot is at index 112:
  2023-08-16 00:47:55.324736  -00deg 48' 46.517584"  -5.176933 arcseconds
  2023-08-16 00:47:55.325735  -00deg 48' 46.516681"  -5.176029 arcseconds
  2023-08-16 00:47:55.326735  -00deg 48' 46.515778"  -5.175126 arcseconds
Worst overshot is at index 70:
  2023-05-26 06:09:37.877960  -00deg 48' 42.952997"   2.818692 arcseconds
  2023-05-26 06:09:37.878960  -00deg 48' 42.952173"   2.819516 arcseconds
  2023-05-26 06:09:37.879960  -00deg 48' 42.951349"   2.820339 arcseconds
'''

from datetime import date, datetime, timedelta
from math import atan, degrees, tau
from skyfield import VERSION, almanac, api, functions
from skyfield.api import load, Topos, wgs84, N, S, E, W
from skyfield.nutationlib import iau2000b

# Avoid 'This probably means that Tcl wasn't installed properly' error.
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

ts = load.timescale()
t0 = ts.utc(2023, 8, 16)
t1 = t0 + 1.0

eph = load('de421.bsp')
earth = eph['earth']
moon = eph['moon']
observer = earth + wgs84.latlon(70, 0)

trace = []
def tracing_at(t):
    trace.append(t)
    return real_at(t)

real_at = observer.at
observer.at = tracing_at

almanac.find_risings(observer, moon, t0, t1)
target = moon

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

times = trace[1:]
for i, t in enumerate(times, 1):
    alt, az, distance = observer.at(t).observe(target).apparent().altaz()
    print(t.utc_strftime(), alt, az)
    ax.plot(az.degrees, alt.degrees, 'o')
    ax.annotate(f'{i}', (az.degrees, alt.degrees), (-1, 1),
                textcoords='offset fontsize')

one_hour = 1.0 / 24.0
one_minute = one_hour / 60
t0 = ts.tt_jd(min(t.tt[0] for t in times) - 10*one_minute)
t1 = ts.tt_jd(max(t.tt[0] for t in times) + 10*one_minute)
t = ts.linspace(t0, t1, 500)
alt, az, distance = observer.at(t).observe(target).apparent().altaz()
ax.plot(az.degrees, alt.degrees, '-')

ax.grid()
fig.savefig('tmp.png')
