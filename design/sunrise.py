from skyfield import almanac, api
from skyfield.api import load, tau

# http://slittlefair.staff.shef.ac.uk/teaching/phy115/session3/page4/page6/page6.html

# sin ALT = sin LAT * sin DEC + cos LAT * cos DEC * cos HA
# sin LAT * sin DEC + cos LAT * cos DEC * cos HA = sin ALT
# cos LAT * cos DEC * cos HA = sin ALT - sin LAT * sin DEC
# cos HA = (sin ALT - sin LAT * sin DEC) / (cos LAT * cos DEC)
# HA = arrcos((sin ALT - sin LAT * sin DEC) / (cos LAT * cos DEC))

#cos HA = (sin ALT - sin LAT sin DEC) / (cos LAT cos DEC)

from numpy import arccos, sin, cos

def _sunrise_hour_angle_radians(latitude, declination, altitude_radians):
    lat = latitude.radians
    dec = declination.radians
    ha = arccos((sin(altitude_radians) - sin(lat) * sin(dec))
                / (cos(lat) * cos(dec)))
    return ha

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
    t = ts.tt_jd(np.arange(start_time.tt, end_time.tt, 0.8))

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

    # We want to return each rising exactly once, so where there are
    # runs of several times `t` that all precede the same rising, let's
    # throw the first few out and keep only the last one.
    i, = np.nonzero(np.diff(difference) < 0.0)

    # When might have rising actually have taken place?  Let's guess by
    # drawing a line between the two (time, hour angles) coordinates and
    # finding where they intersect zero.
    a = tau - difference[i]
    b = difference[i + 1]
    tt = t.tt
    new_tt = (b * tt[i] + a * tt[i+1]) / (a + b)
    t = ts.tt_jd(new_tt)

    ha_per_day = tau

    for i in 0, 1: #2, 3:
        fastify(t)
        ha, dec, _ = observer.at(t).observe(target).apparent().hadec()
        desired_ha = - _sunrise_hour_angle_radians(latitude, dec, h)
        ha_adjustment = desired_ha - ha.radians
        print('ha_adjustment:', ha_adjustment)
        timebump = ha_adjustment / ha_per_day
        print('Timebump(s):', timebump * 24.0 * 3600.0)
        t = ts.tt_jd(t.whole, t.tt_fraction + timebump)

    # Test of how we did.
    t_final = t
    alt, az, _ = observer.at(t_final).observe(target).apparent().altaz()
    diff = horizon_degrees - alt.degrees
    print('New alts:', diff.min(), 'to', diff.max(), 'degrees')
    print('Detailed differences (degrees):', diff)

    return t

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

    if 0:
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
