#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import date, datetime, timedelta
from math import atan, degrees, tau
from skyfield import VERSION, almanac, api, functions
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

latlon = wgs84.latlon(latitude, 0.0 * E)
observer = earth + latlon
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
arrow_style = {
    'arrowstyle': 'simple',
    'fc': 'black', 'lw': 0.2, 'relpos': (0, 0.5),
    'shrinkA': 4, 'shrinkB': 12,
}

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
    arrow_x = 30.0
    arrow_y = 10.0
    for i, t in enumerate(traced_times):
        print(t.utc_strftime())
        if i == 0:
            continue
        alt, az, distance = observer.at(t).observe(target).apparent().altaz()
        ax.plot(az.degrees, alt.degrees, 'o')
        arrow_y -= 9.0
        print((arrow_x, arrow_y))
        ax.annotate(
            f'{i}  ({t[0].utc_strftime("%H:%M:%S")})',
            xy=(az.degrees, alt.degrees), xytext=(arrow_x, arrow_y),
            textcoords='offset points',
            arrowprops=arrow_style,
        )

    if 0: #target is moon:
        # Trying to find the intersection by hand:
        ty = traced_times[-2][0]
        tz = traced_times[-1][0]
        alty, azy, distance = observer.at(ty).observe(target).apparent().altaz()
        altz, azz, distance = observer.at(tz).observe(target).apparent().altaz()

        ax.plot([azy.degrees, azz.degrees],
                [alty.degrees, altz.degrees], '--r')

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
        #arrow_angle -= 20.0 * degree
        # ax.annotate(str('end'), (az.degrees, alt.degrees), (arrow_x, arrow_y),
        #             'data', 'offset points', arrow_style)

    ty = traced_times[-2][0]
    tz = traced_times[-1][0]

    if 1: #target is moon:  # try number two
        # Trying to find the intersection by hand:
        vy = observer.at(ty).observe(target).apparent().frame_xyz(latlon).km
        vz = observer.at(tz).observe(target).apparent().frame_xyz(latlon).km
        print(vy)
        print(vz)
        vy /= functions.length_of(vy)
        vz /= functions.length_of(vz)
        print(vy)
        print(vz)
        print('target alt:', -0.0147393)
        x,y,z = vy
        print('y alt:', np.atan2(z, np.sqrt(x*x + y*y)))
        x,y,z = vz
        print('z alt:', np.atan2(z, np.sqrt(x*x + y*y)))

        # so let's try drawing the line between those two points
        px = []
        py = []
        p_vector = np.linspace(0, 1, 100)
        for p in p_vector:
            vp = vy * (1 - p) + vz * p
            x,y,z = vp
            p_alt = np.atan2(z, np.sqrt(x*x + y*y)) / tau * 360.0
            p_az = np.atan2(y, x) / tau * 360.0
            print(vp, p_alt, p_az)
            px.append(p_az)
            py.append(p_alt)

        ax.plot(px, py, '--g')

        #print('So which point on the line between those is at target alt?')

        # BUT: how would that help? I need to know time, and the parameter
        # `t` of the great-circle-segment swept out by the formula does NOT
        # move across the sphere at a uniform rate!  So even if I answered
        # the above question, I wouldn't know the corresponding time.

    # Drat!  The great-circle line drawn above misses the rising point
    # on the horizon, just like the naive line drawn directly across the
    # plot.  So the curving path of the Moon isn't just an artifact of
    # the projection.  Is the Moon's ra/dec velocity noticably different
    # from one end of the line to the other?  And how big is its ra/dec
    # velocity compared to its velocity in the sky overall?

    # A: Because the Moon isn't at 0° declination, on the equator, so
    # the path it's taking isn't a great circle; it's a little one; and
    # even on a small local scale, a little circle has a tighter curve
    # than a great circle.

    if target is moon:
        py = observer.at(ty).observe(target).apparent()
        pz = observer.at(tz).observe(target).apparent()
        pyk = py.velocity.km_per_s
        pzk = pz.velocity.km_per_s
        print(pyk)
        print(pzk)
        print((pyk - pzk) / functions.length_of(pyk))

    ax.set(xlabel='Azimuth (°)', ylabel='Altitude (°)',
           title='The search for the moment of rising on 2023-02-19')
    ax.axhline(horizon)
    ax.axvline(180.0)
    return fig

if 1:
    #fig, (ax, ax2) = plt.subplots(2)
    fig, ax = plt.subplots()
    t0 = ts.utc(2023, 2, 19)
    t = ts.linspace(
        ts.utc(2023, 2, 19, 11, 21),
        #ts.utc(2023, 2, 19, 11, 5),
        ts.utc(2023, 2, 19, 11, 36),
        299,
    )
    #alt, az, distance = observer.at(t).observe(moon).apparent().altaz()
    ha, dec, distance = observer.at(t).observe(moon).apparent().hadec()
    minute = (t - t0) * 1440 - 11*60
    ax.plot(minute, ha._degrees)

    altitude_radians = -0.01473935
    rha = almanac._rising_hour_angle(latlon.latitude, dec, altitude_radians)
    ax.plot(minute, rha / tau * 360.0)
    ax.set(xlabel='Time (minutes into the hour)',
           ylabel='Hour Angle (°)')

    ax.set(xlabel='Time (minutes into the hour)',
           ylabel='Hour Angle (°)')

    # ax2.plot(minute, dec.degrees)
    fig.savefig('tmp4.png')

# make_plot(stars[2], horizon).savefig('tmp1.png')
# make_plot(stars[1], horizon).savefig('tmp2.png')
make_plot(moon).savefig('tmp3.png')
