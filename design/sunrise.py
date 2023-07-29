from skyfield import almanac, api
from skyfield.api import load, tau

# sunrise_t, sunset_t = t

# app = (earth + bluffton).at(sunrise_t).observe(sun).apparent()
# ha, dec, distance = app.hadec()
# alt, az, distance = app.altaz()

# print('Dec in degrees at sunrise:', dec.degrees)
# print('Alt in degrees at sunrise:', alt.degrees)
# print('HA in hours at sunrise:', ha.hours)

# def f(t, msg):
#     app = (earth + bluffton).at(t).observe(sun).apparent()
#     ha, dec, distance = app.hadec()
#     alt, az, distance = app.altaz()
#     print(f'HA in hours {msg}:', ha.hours)

# f(sunset_t, 'at sunset')
# f(ts.tt_jd(sunrise_t.tt - 1/24.0/3600.0), 'at sunrise - 1 second')
# f(ts.tt_jd(sunrise_t.tt + 1/24.0/3600.0), 'at sunrise + 1 second')

# http://slittlefair.staff.shef.ac.uk/teaching/phy115/session3/page4/page6/page6.html

# sin ALT = sin LAT * sin DEC + cos LAT * cos DEC * cos HA
# sin LAT * sin DEC + cos LAT * cos DEC * cos HA = sin ALT
# cos LAT * cos DEC * cos HA = sin ALT - sin LAT * sin DEC
# cos HA = (sin ALT - sin LAT * sin DEC) / (cos LAT * cos DEC)
# HA = arrcos((sin ALT - sin LAT * sin DEC) / (cos LAT * cos DEC))

#cos HA = (sin ALT - sin LAT sin DEC) / (cos LAT cos DEC)

from numpy import arccos, sin, cos

# alt, az, _ = app.altaz()

# ALT =
# print('ALT', alt.radians)
# DEC = dec.radians
# print('DEC', dec.radians)
# LAT = bluffton.latitude.radians
# print('LAT', bluffton.latitude.radians)

def _sunrise_hour_angle_radians(latitude, declination, altitude_radians):
    lat = latitude.radians
    dec = declination.radians
    ha = arccos((sin(altitude_radians) - sin(lat) * sin(dec))
                / (cos(lat) * cos(dec)))
    return ha

# ha = _sunrise_hour_angle_radians(bluffton.latitude, dec, stdalt)

# print('HA, computed:', ha / tau * 24.0)

import numpy as np

from skyfield.nutationlib import iau2000b_radians
def fastify(t):
    t._nutation_angles_radians = iau2000b_radians(t)

def find_sunrise(observer, target, horizon_degrees, start_time, end_time):
    geo = observer.vector_functions[-1]
    latitude = geo.latitude
    h = horizon_degrees / 360.0 * tau

    # Build an array of times 0.8 days apart from start_time to end_time.
    ts = start_time.ts
    tt0 = start_time.tt
    tt1 = end_time.tt
    t = ts.tt_jd(np.arange(tt0, tt1, 0.8))

    # Determine the target's hour angle and declination at those times.
    fastify(t)
    ha, dec, _ = observer.at(t).observe(target).apparent().hadec()

    # Invoke our geometry formula: for each time `t`, predict the hour
    # angle at which the target will next reach the horizon, if its
    # declination were to remain constant.
    setting_ha = _sunrise_hour_angle_radians(latitude, dec, h)
    rising_radians = - setting_ha

    # So at each time `t`, how many radians is the target's hour angle
    # from the target's next 'ideal' rising?
    difference = ha.radians - rising_radians
    difference %= tau
    print('HA jump (deg):', difference / tau * 360.0)

    # We want to return each rising exactly once, so where there are
    # runs of several times `t` that all precede the same rising, let's
    # throw the first few out and keep only the last one.
    i, = np.nonzero(np.diff(difference) < 0.0)
    print(i)

    # When might have rising actually have taken place?  Let's guess by
    # drawing a line between the two (time, hour angles) coordinates and
    # finding where they intersect zero.
    a = tau - difference[i]
    b = difference[i + 1]
    tt = t.tt
    new_tt = (b * tt[i] + a * tt[i+1]) / (a + b)

    # Time to check how good those guesses were!  Let's plug those TT
    # times in and see where exactly the target really was.
    t2 = ts.tt_jd(new_tt)
    fastify(t2)
    ha2, dec2, _ = observer.at(t2).observe(target).apparent().hadec()

    # (What is the diff between expected dec and actual dec?  Is that
    # what makes all the difference?)
    print('Dec diff (deg):', dec2.degrees - dec.degrees[i])

    # The target won't have been exactly at the right place because its
    # declination won't really have stayed constant, so let's repeat the
    # above computation to bring us one step closer.
    setting_ha2 = _sunrise_hour_angle_radians(latitude, dec2, h)
    rising_ha2 = - setting_ha2

    adjustment = rising_ha2 - ha2.radians  # radians
    rise = ha2.radians - ha.radians[i]     # radians
    run = t2.tt - t[i].tt                  # TT
    slope = rise / run                     # radians actually turned per TT
    timebump = adjustment / slope          # rad * (rad / TT) -> TT

    print('Timebump(s):', timebump * 24.0 * 3600.0)

    t3 = ts.tt_jd(t2.tt + timebump)

    fastify(t3)
    ha3, dec3, _ = observer.at(t3).observe(target).apparent().hadec()
    setting_ha3 = _sunrise_hour_angle_radians(latitude, dec3, h)
    rising_ha3 = - setting_ha3

    adjustment = rising_ha3 - ha3.radians  # radians
    rise = ha3.radians - ha2.radians       # radians
    run = t3.tt - t2.tt                    # TT
    slope = rise / run                     # radians actually turned per TT
    timebump = adjustment / slope          # rad * (rad / TT) -> TT

    t4 = ts.tt_jd(t3.tt + timebump)

    # Test of how we did.
    t_final = t4
    alt, az, _ = observer.at(t4).observe(target).apparent().altaz()
    diff = horizon_degrees - alt.degrees
    print('New alts:', diff.min(), 'to', diff.max(), 'degrees')
    print('Detailed differences (degrees):', diff)

    return t_final

    # adjustment = stdalt - alt.radians
    # print(max(abs(adjustment)) / tau * 24.0 * 3600.0,
    #       'radians adj, in seconds of day (approx)')
    # print(max(adjustment / tau * 24.0 * 3600.0),
    #       'max')

    #ha, dec, _ = observer.at(t).observe(target).apparent().hadec()

#find_sunrise(earth + bluffton, sun, stdalt, t0, t1)

def main():
    ts = load.timescale()
    eph = load('de421.bsp')
    earth = eph['earth']
    #bluffton = earth + api.wgs84.latlon(40.8939, -83.8917)
    bluffton = api.wgs84.latlon(40.8939, -83.8917)

    t0 = ts.utc(2020, 1, 1)
    t1 = ts.utc(2020, 1, 10)
    #t1 = ts.utc(2021, 1, 1)
    # t0 = ts.utc(2018, 9, 12, 4)
    # t1 = ts.utc(2018, 9, 13, 4)

    from time import time

    if 1:
        sun = eph['sun']
        stdalt = -0.8333 #/ 360.0 * tau

        T0 = time()
        f = almanac.sunrise_sunset(eph, bluffton)
        t_old, y = almanac.find_discrete(t0, t1, f)
        DUR_OLD = time() - T0

    else:
        sun = eph['mercury']
        stdalt = -0.8333 #/ 360.0 * tau

        T0 = time()
        f = almanac.risings_and_settings(
            eph, sun, bluffton,
            horizon_degrees=-stdalt,
            radius_degrees=0.0,
        )
        #f = almanac.sunrise_sunset(eph, bluffton)
        t_old, y = almanac.find_discrete(t0, t1, f)
        DUR_OLD = time() - T0

    print('How good are old find_discrete results?')
    alt, az, _ = (earth + bluffton).at(t_old).observe(sun).apparent().altaz()
    diff_degrees = alt.degrees - stdalt
    print(diff_degrees)
    print()

    # print(t.utc_iso())
    # print(y)

    T0 = time()
    t_new = find_sunrise(earth + bluffton, sun, stdalt, t0, t1)
    DUR_NEW = time() - T0

    t_old = t_old[::2]

    print(len(t_old), 'old')
    print(len(t_new), 'new')

    diff_seconds = (t_old.tt - t_new.tt) * 24.0 * 3600.0
    print('Max difference (seconds):', max(
        abs(diff_seconds)
    ))

    alt, az, _ = (earth + bluffton).at(t_new).observe(sun).apparent().altaz()
    diff_degrees = alt.degrees - stdalt

    if 1:
        import matplotlib.pyplot as plt
        fig, (ax, ax2) = plt.subplots(2)
        ax.plot(t_new.J, diff_seconds, label='diff (sec)', linestyle='--')
        ax.grid()
        ax2.plot(t_new.J, diff_degrees, label='altitude °', linestyle='--')
        ax2.grid()
        #ax.set(xlabel='time (s)', ylabel='voltage (mV)', title='Title')
        #ax.set_aspect(aspect=1.0)
        #ax.axhline(1.0)
        #ax.axvline(1.0)
        #plt.legend()
        fig.savefig('tmp.png')

    print(DUR_OLD, DUR_NEW, DUR_OLD / DUR_NEW)

    #print('HA:', HA / tau * 24.0 - 12)
    #print(HA - api.tau)

main()
