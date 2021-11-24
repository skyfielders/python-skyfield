from skyfield import almanac, api
from skyfield.api import load, tau

ts = load.timescale()
eph = load('de421.bsp')
sun = eph['sun']
earth = eph['earth']
#bluffton = earth + api.wgs84.latlon(40.8939, -83.8917)
bluffton = api.wgs84.latlon(40.8939, -83.8917)

t0 = ts.utc(2020, 1, 1)
t1 = ts.utc(2021, 1, 1)
# t0 = ts.utc(2018, 9, 12, 4)
# t1 = ts.utc(2018, 9, 13, 4)

# print(t.utc_iso())
# print(y)

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

stdalt = -0.833 / 360.0 * tau
# ha = _sunrise_hour_angle_radians(bluffton.latitude, dec, stdalt)

# print('HA, computed:', ha / tau * 24.0)

import numpy as np

def find_sunrise(observer, body, horizon, start_time, end_time):
    ts = start_time.ts
    tt0 = start_time.tt
    tt1 = end_time.tt
    t = ts.tt_jd(np.arange(tt0, tt1, 0.8))
    ha, dec, _ = observer.at(t).observe(body).apparent().hadec()
    def print(*args):
        pass
    print(ha.hours)
    print(dec.degrees)

    geo = observer.vector_functions[-1]
    setting_ha = _sunrise_hour_angle_radians(geo.latitude, dec, stdalt)

    rising_radians = - setting_ha
    difference = ha.radians - rising_radians
    print(difference)
    # difference = np.unwrap(difference)
    difference %= tau
    print(difference)

    i, = np.nonzero(np.diff(difference) < 0.0)
    print(i)
    a = tau - difference[i]
    b = difference[i + 1]
    print(a)
    print(b)
    tt = t.tt
    new_tt = (b * tt[i] + a * tt[i+1]) / (a + b)
    print(tt)
    print(new_tt)

    t2 = ts.tt_jd(new_tt)
    ha2, dec2, _ = observer.at(t2).observe(body).apparent().hadec()
    setting_ha2 = _sunrise_hour_angle_radians(geo.latitude, dec2, stdalt)
    rising_ha2 = - setting_ha2
    # alt, az, _ = observer.at(t).observe(body).apparent().altaz()
    # print(alt.degrees, 'alt')
    # print(alt.radians, 'alt radians')
    # print(stdalt, 'desired radians')

    # try a->t ?

    adjustment = rising_ha2 - ha2.radians
    print(max(abs(adjustment)) / tau * 24.0 * 3600.0,
          'max adjustment (~seconds of day)')

    rise = ha2.radians - ha.radians[i]
    run = t2.tt - t[i].tt
    slope = rise / run

    timebump = adjustment / slope
    print(timebump)

    t3 = ts.tt_jd(t2.tt + timebump)
    ha3, dec3, _ = observer.at(t3).observe(body).apparent().hadec()
    setting_ha3 = _sunrise_hour_angle_radians(geo.latitude, dec3, stdalt)
    rising_ha3 = - setting_ha3

    print(ha3.radians - rising_ha3)
    print(max(abs((ha3.radians - rising_ha3))) / tau * 24.0 * 3600.0,
          'seconds')

    return t3

    # adjustment = stdalt - alt.radians
    # print(max(abs(adjustment)) / tau * 24.0 * 3600.0,
    #       'radians adj, in seconds of day (approx)')
    # print(max(adjustment / tau * 24.0 * 3600.0),
    #       'max')

    #ha, dec, _ = observer.at(t).observe(body).apparent().hadec()

#find_sunrise(earth + bluffton, sun, stdalt, t0, t1)

from time import time
T0 = time()
t_new = find_sunrise(earth + bluffton, sun, stdalt, t0, t1)
DUR_NEW = time() - T0

T0 = time()
t_old, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(eph, bluffton))
DUR_OLD = time() - T0

t_old = t_old[::2]

print(len(t_old), 'old')
print(len(t_new), 'new')

print(max(abs(t_old.tt - t_new.tt)))
print(DUR_OLD, DUR_NEW, DUR_OLD / DUR_NEW)

#print('HA:', HA / tau * 24.0 - 12)
#print(HA - api.tau)
