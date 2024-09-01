#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import date, datetime, timedelta
from math import atan, degrees
from skyfield import VERSION, almanac, api
from skyfield.api import load, Topos, wgs84, N, S, E, W
from skyfield.nutationlib import iau2000b

import numpy as np
import matplotlib.pyplot as plt

latitude = 70.0 * N

eph = load('de421.bsp')
earth   = eph['earth']
moon    = eph['moon']

ts = load.timescale()
t0 = ts.utc(2023, 2, 19)
t1 = t0 + 1.0

#['2023-02-19 11:35:05 UTC'] vs 11:27:58Z

trace = []

def patched_at(t):
    trace.append(t)
    return real_observer_at(t)

observer = earth + wgs84.latlon(latitude, 0.0 * E)
real_observer_at = observer.at
observer.at = patched_at

ra, dec, distance = observer.at(t0).observe(moon).apparent().radec('date')
print(ra.hours, dec.degrees)
ra, dec, distance = observer.at(t1).observe(moon).apparent().radec('date')
print(ra.hours, dec.degrees)

t = ts.utc(2023, 2, 19, 11, range(15, 55))

star_declinations = [-21.0, -20.95, -20.9, -20.85, -20.8]
stars = [api.Star(ra_hours=21.5, dec_degrees=d) for d in star_declinations]

star_curves = []
for star in stars:
    alt, az, distance = observer.at(t).observe(star).apparent().altaz()
    star_x = az.degrees
    star_y = alt.degrees
    star_curves.append([star, star_x, star_y])

alt, az, distance = observer.at(t).observe(moon).apparent().altaz()
moon_x = az.degrees
moon_y = alt.degrees

horizon = -0.8444219217855957  # later, compute this instead of hard-coding
arrow_style = {'arrowstyle': 'simple', 'lw': 0.2, 'fc': 'black'}

def make_plot(target, horizon_degrees=None):
    trace.clear()
    t, y = almanac.find_risings(observer, target, t0, t1, horizon_degrees)
    traced_times = trace[:]

    fig, ax = plt.subplots()

    for star, star_x, star_y in star_curves:
        ax.plot(star_x, star_y, 'b')
        ax.annotate(f'{star.dec.degrees}°', (star_x[-1], star_y[-1]), (-23, -9),
                    'data', 'offset points')

    if target is moon:
        ax.plot(moon_x, moon_y, '.')

    ax.grid()

    degree = api.tau / 360.0
    arrow_angle = 80.0 * degree
    arrow_length = 25.0
    for i, t in enumerate(traced_times):
        print(t.utc_strftime())
        if i == 0:
            continue
        alt, az, distance = observer.at(t).observe(target).apparent().altaz()
        ax.plot(az.degrees, alt.degrees, 'o')
        arrow_angle -= 20.0 * degree
        arrow_x = np.cos(arrow_angle) * arrow_length
        arrow_y = np.sin(arrow_angle) * arrow_length
        ax.annotate(f'{i}  ({t[0].utc_strftime("%H:%M:%S")})',
                    (az.degrees, alt.degrees), (arrow_x, arrow_y),
                    'data', 'offset points', arrow_style)

    if 0: #target is moon:
        # Trying to find the intersection by hand:
        ty = traced_times[-2][0]
        tz = traced_times[-1][0]
        alty, az, distance = observer.at(ty).observe(target).apparent().altaz()
        altz, az, distance = observer.at(tz).observe(target).apparent().altaz()

        ay = alty.radians - -0.0147393
        az = altz.radians - -0.0147393
        t = ty - ay / ((az - ay) / (tz - ty))
        # t = ty - ay / (az - ay) * (tz - ty)
        print(ty.utc_strftime(), ay)
        print(tz.utc_strftime(), az)
        print(t.utc_strftime())

        alt, az, distance = observer.at(t).observe(target).apparent().altaz()
        print(alt.degrees, az.degrees)

        ax.plot(az.degrees, alt.degrees, 'o')
        arrow_angle -= 20.0 * degree
        ax.annotate(str('end'), (az.degrees, alt.degrees), (arrow_x, arrow_y),
                    'data', 'offset points', arrow_style)

    ax.set(xlabel='Azimuth (°)', ylabel='Altitude (°)',
           title='The search for the moment of rising on 2023-02-19')
    ax.axhline(horizon)
    ax.axvline(180.0)
    return fig

make_plot(stars[2], horizon).savefig('tmp1.png')
make_plot(stars[1], horizon).savefig('tmp2.png')
make_plot(moon).savefig('tmp3.png')
