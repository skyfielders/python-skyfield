from skyfield import almanac
from skyfield import api

ts = api.load.timescale()
t0 = ts.tt(-1000, 1, 1)
t1 = ts.tt(2000, 1, 1)
days = int(t1 - t0)

if 1:
    t = ts.tt(-1000, 1, range(days))
    gy, gm, gd = t.tt_calendar()[:3]
    ts.julian_calendar_cutoff = 99999999999999999 #api.GREGORIAN_START
    jy, jm, jd = t.tt_calendar()[:3]
    print(gd - jd)

    # import numpy as np
    # t = np.arange(0.0, 2.0, 0.01)
    # s = 1 + np.sin(2 * np.pi * t)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(t.J, gd - jd, '.') #, label='label', linestyle='--')
    ax.axvline(325, color='blue')
    #ax.plot(325, 0, '.') #, label='label', linestyle='--')
    # ax.set(xlabel='time (s)', ylabel='voltage (mV)', title='Title')
    # ax.set_aspect(aspect=1.0)
    ax.grid()
    # plt.legend()
    ax.set_ylim(-10, 15)
    fig.savefig('tmp.png')
    exit()

eph = api.load('de422.bsp')
t, y = almanac.find_discrete(t0, t1, almanac.seasons(eph))

vernal = (y == 0)
t = t[vernal]
y = y[vernal]

#cal = t.ut1_calendar
cal = lambda: t.utc
gday = cal()[2]
for yi, ti in zip(y, t):
    if 325 < ti.J <= 326:
        print(yi, almanac.SEASON_EVENTS[yi], ti.utc_iso(' '))

ts.julian_calendar_cutoff = api.GREGORIAN_START

del t.utc
jday = cal()[2]
for yi, ti in zip(y, t):
    if 325 < ti.J <= 326:
        print(yi, almanac.SEASON_EVENTS[yi], ti.utc_iso(' '))

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(t.J, jday, 'o')
ax.plot(t.J, gday, '.') #, alpha=0.25)
ax.grid()
ax.axvline(325, color='r')
fig.savefig('tmp.png')
