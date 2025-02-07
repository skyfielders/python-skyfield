#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import date, datetime, timedelta
from math import atan, degrees
from skyfield import VERSION, almanac
from skyfield.api import Loader, Topos, wgs84, N, S, E, W
from skyfield.nutationlib import iau2000b

def getHorizon(t):
    # calculate the angle of the moon below the horizon at moonrise/set

    position = earth.at(t).observe(moon)   # at noontime (for daily average distance)
    distance = position.apparent().radec(epoch='date')[2]
    dist_km = distance.km
    sdm = degrees(atan(1737.4/dist_km))   # volumetric mean radius of moon = 1737.4 km
    horizon = sdm + 0.5666667	# moon's equatorial radius + 34' (atmospheric refraction)
    return horizon

def f_moon(topos, degBelowHorizon):
    # Build a function of time that returns the moon state.
    topos_at = (earth + topos).at

    def is_moon_up_at(t):
        """The function that this returns will expect a single argument that is a 
		:class:`~skyfield.timelib.Time` and will return ``True`` if the moon is up,
		else ``False``."""
        t._nutation_angles = iau2000b(t.tt)
        # Return `True` if the moon has risen by time `t`.
        return topos_at(t).observe(moon).apparent().altaz()[0].degrees > -degBelowHorizon

    is_moon_up_at.rough_period = 0.5  # twice a day
    return is_moon_up_at

d = date(2023, 2, 18)       # modify the starting date as required
lat = 70.0                  # modify latitude as required (use float for Topos)

load = Loader("./")
ts = load.timescale()
eph = load('de421.bsp')
earth   = eph['earth']
moon    = eph['moon']
dt = datetime(d.year, d.month, d.day, 0, 0, 0)
print("Skyfield version = {};  Results for latitude {} ...".format(VERSION,lat))

for i in []: #range(3):      # sample 3 successive dates

    t0     = ts.ut1(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    t0noon = ts.ut1(dt.year, dt.month, dt.day, dt.hour+12, dt.minute, dt.second)
    t1     = ts.ut1(dt.year, dt.month, dt.day+1, dt.hour, dt.minute, dt.second)

    topos = wgs84.latlon(lat, 0.0 * E, elevation_m=0.0)
    observer = earth + topos
    horizon = getHorizon(t0noon)

    moonrise, yR = almanac.find_risings(observer, moon, t0, t1, -horizon)  # -horizon-0.048 fixes problem
    moonset,  yS = almanac.find_settings(observer, moon, t0, t1, -horizon)

    locn = Topos(lat, "0.0 E")
    moonevent, y = almanac.find_discrete(t0, t1, f_moon(locn, horizon))

    print("{}:  {} x moonrise, {} x moonset;  {} x moonevent;  horizon = {:.4f}°".format(dt.strftime("%Y-%m-%d"),len(yR),len(yS),len(y),horizon))
    print("   moonrise ",moonrise.utc_iso(' '),"   ",yR)
    print("   moonset  ",moonset.utc_iso(' ') ,"   ",yS)
    print("   moonevent",moonevent.utc_iso(' '),"   ",y)
    
    dt += timedelta(days=1)

topos = wgs84.latlon(lat, 0.0 * E, elevation_m=0.0)
observer = earth + topos
t0 = ts.utc(2023, 2, 19)
t1 = ts.utc(2023, 2, 20)

horizon = - getHorizon(t0)
print('Horizon:', horizon)

t, y = almanac.find_risings(observer, moon, t0, t1, horizon)
print(t.utc_strftime(), 'vs 11:27:58Z')
print(y)

import numpy as np
t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2 * np.pi * t)

t = ts.linspace(t0, t1, 3600)
a = observer.at(t).observe(moon).apparent()
alt, az, distance = a.altaz()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, [ax, ax2] = plt.subplots(2)
#ax.plot(t.tt - t0.tt, alt.degrees, label='Altitude') #, linestyle='--')
ax.plot(az.degrees, alt.degrees, label='Altitude') #, linestyle='--')
ax.grid()
#ax.set(xlabel='time (s)', ylabel='voltage (mV)', title='Title')
#ax.set_aspect(aspect=1.0)
ax.axhline(horizon)
ax.axvline(180.0)
#; ax.axvline(1.0); ax.set_xlim(10, 20)
#ax.set_ylim(-.850, -.8)
#ax.set_ylim(-.850, -.75)
ax.set_ylim(-1, -.75)
ax.set_xlim(170.0, 190.0)

ra, dec, distance = a.radec('date')
ax2.plot(az.degrees, dec.degrees)
ax2.set_xlim(170.0, 190.0)
ax2.set_ylim(-21, -20)
ax2.grid()

#plt.legend()

cutoff = None  # -2
t2 = ts.tt_jd(almanac.TRACE[:cutoff])
alt, az, distance = observer.at(t2).observe(moon).apparent().altaz()

ax.plot(az.degrees, alt.degrees, '.', label='Tries')

fig.savefig('tmp.png')

