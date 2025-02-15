#!/usr/bin/env python3
"""Plot a diagram of a difficult search for Moonrise.

The `design/test_sunrise_moonrise.py` finishes with a search for high
latitude (70°N) moonrises, and prints out the worst miss with a negative
elevation versus the horizon, and the worst positive:

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

Each line is one millisecond apart.  In each case we really wanted the
search to find the millisecond where the arcseconds flip from negative
to positive, but missed.

Why?

This script takes the 2023-08-16 case shown above, instruments the .at()
function to record the calls made by the search for moonrise, and plots
the result so we can see what's happening.  I'm planning to keep this
script around in the repository, so that if any future users find fault
with the new risings routine, I can quickly adapt this script to their
case to see what's going on.

"""
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

# Tweak these.

TARGET = 'moon'
YEAR_MONTH_DAY = 2023, 8, 16
ZOOM_ARCSECONDS = None

# You will get a diagram in return.

ts = load.timescale()
t0 = ts.utc(*YEAR_MONTH_DAY) #- 0.1
t1 = t0 + 1.0 #+ 1.1

eph = load('de421.bsp')
earth = eph['earth']
target = eph[TARGET]
observer = earth + wgs84.latlon(70, 0)

trace = []
def tracing_at(t):
    trace.append(t)
    return real_at(t)

real_at = observer.at
observer.at = tracing_at

final_t, y = almanac.find_risings(observer, target, t0, t1)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

# Drop trace[0], since it's different: it's the initial big-picture
# search for intervals that have risings in them at all, which are then
# filtered out and refined through the subsequent searches.
times = trace[1:]

# And add the final answer back in as though it were a final iteration.
times.append(final_t)

for i, t in enumerate(times, 1):
    alt, az, distance = observer.at(t).observe(target).apparent().altaz()
    print(t.utc_strftime('%Y-%m-%d %H:%M:%S.%f'), alt, az)
    ax.plot(az.degrees, alt.degrees, 'o')
    ax.annotate(f'{i}', (az.degrees, alt.degrees), (-1 + 2 * (i%2), 1),
                textcoords='offset fontsize')

h = almanac.build_horizon_function(target)
horizon_degrees = h(distance) / tau * 360.0
horizon_degrees = horizon_degrees[0]
print(f'Horizon: {horizon_degrees}°')
ax.axhline(horizon_degrees)

miss_degrees = alt.degrees[0] - horizon_degrees
miss_arcseconds = miss_degrees * 3600.0
print(f'Error: {miss_arcseconds} arcseconds')

if ZOOM_ARCSECONDS is not None:
    ax.set_ylim(
        alt.degrees[0] - ZOOM_ARCSECONDS / 3600.0,
        alt.degrees[0] + ZOOM_ARCSECONDS / 3600.0,
    )

one_hour = 1.0 / 24.0
one_minute = one_hour / 60
t0 = ts.tt_jd(min(t.tt[0] for t in times) - 10*one_minute)
t1 = ts.tt_jd(max(t.tt[0] for t in times) + 10*one_minute)
t = ts.linspace(t0, t1, 500)
alt, az, distance = observer.at(t).observe(target).apparent().altaz()
ax.plot(az.degrees, alt.degrees, '-')

ax.grid()
fig.savefig('tmp.png')
print('Wrote tmp.png')
