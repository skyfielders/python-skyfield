from skyfield import api
from skyfield.api import load, T0
from skyfield.constants import tau
from skyfield.data.mpc import *

time = __import__('time').time

with load.open(COMET_URL) as f:
    print(f)
    t0 = time()
    df = load_comets_dataframe_slow(f)
    print(time() - t0)

print(df.head())

# t0 = time()
# d = df['designation']
# e = d.str[-7]
# print(time() - t0)
# #print(e)

# for dd in df['designation'].values:
#     if dd.count('K') > 1:
#         print(dd)

# print(df.ix[0])
# i = df.info()
# print(i)
# i = df.ix[0]
# print(i)
# # ref?
# #print(df.reference.value_counts())
# # 'computer_name',
# print(time() - t0)

from skyfield.keplerlib import KeplerOrbit

ts = load.timescale(builtin=True)

# df = df[df.designation == 1]

from numpy import sqrt

print('Test:',
      186, 'vs',
      0.914 / (1.0 - 0.995086))

df['semimajor_axis_au'] = (
    df['perihelion_distance_au'] / (1.0 - df['eccentricity'])
)

# Periapsis distance = a(1-e)
# Periapsis distance / (1-e) = a
# Periapsis distance / (1-e) = semi-major axis
# semi-major axis = Periapsis distance / (1-e)

# df['mean_anomaly_degrees'] = 0.40465

# comet = df.iloc[0:1]
# print(comet.ix[0])

from skyfield.data.gravitational_parameters import GM_dict
from skyfield.constants import AU_KM, DAY_S

mu_km3_s2 = GM_dict[10]
print('mu_km3_s2:', mu_km3_s2)  # [m3 s−2]  132712440041.93938  11 + 9 = 20  OK
mu_au3_d2 = mu_km3_s2 / (AU_KM**3.0) * (DAY_S**2.0)  # ok?
print('mu_au3_d2:', mu_au3_d2)

row = df.ix[0]

t_perihelion = ts.tt(
    row.perihelion_year, row.perihelion_month, row.perihelion_day
)

df['mean_anomaly_degrees'] = (
    sqrt(mu_au3_d2 / (row.semimajor_axis_au ** 3.0))
    *
    (ts.J2000.tt - t_perihelion.tt)
    * 360.0 / tau
)

print('Days:', ts.J2000.tt - t_perihelion.tt)
print('Back-of-envelope M:', 1007 / (2532 * 365.2425) * 360.0)
print('Back-of-envelope M:', 1007 / (2533 * 365.2425) * 360.0)
print('Back-of-envelope M:', 1007 / (2534 * 365.2425) * 360.0)

print(
    '================== Semimajor:', row.semimajor_axis_au,
)

print(
    '================== Candidate orbital period (years)?:',
    360.0 / sqrt(mu_au3_d2 / (row.semimajor_axis_au ** 3.0))
    / 365.2425 / tau
)

print(
    '================== Candidate M mean_anomaly:',
    sqrt(mu_au3_d2 / (row.semimajor_axis_au ** 3.0))
    *
    (ts.J2000.tt - t_perihelion.tt)
    * 360.0 / tau #??
)

comet = df.iloc[0:1]

k = KeplerOrbit.from_dataframe(comet, ts)

from skyfield.data.spice import inertial_frames

# t = ts.tt_jd(2450538.4378482755)
# print(t.utc_jpl())

print(comet.ix[0]['designation'])

print('CometEls.txt gives perihelion as: 1997 03 29.6333')
print('HORIZONS gives perihelion as:', ts.tt_jd(2450538.4378482755).utc_jpl())

from math import sqrt
#print(sqrt(3.552056259323088E+00**2 + 8.386813198052002E-01**2 + 4.323481144417487E+01**2))
# print('Cut and paste from HORIZONS distance:', sqrt(1.777310651689592E+00**2 + 1.638390146876578E+00**2 + 2.712743223120575E+01**2))

# for month in range(3, 8):
#     t = ts.tdb(2008, month - 9, 15)
#     #t = ts.tdb(2020, month, 24)
#     print('====== Today ======', t.utc_jpl())
#     p = k.at(t)
#     ra, dec, distance = p.radec()

#     from skyfield.functions import length_of
#     print(p.position.au)
#     print('Distance from Sun (AU):', length_of(p.position.au))

eph = api.load('de421.bsp')
t = ts.utc(2020, 5, 31)
#t = ts.utc(2020, 6, 2)

if True:
    k._rotation = inertial_frames['ECLIPJ2000'].T
else:
    # (Reduces accuracy:)
    from skyfield.functions import _mxm, rot_x
    epoch = t
    _, d_eps = epoch._nutation_angles_radians
    true_obliquity = epoch._mean_obliquity_radians + d_eps
    k._rotation = _mxm(rot_x(- true_obliquity), epoch.M).T  # mxmxv?
    #position_au = _mxv(rotation, position_au)

p = eph['earth'].at(t).observe(k)#.apparent()
p.position.au += eph['sun'].at(t).position.au
ra, dec, distance = p.radec()

print(t.utc_iso(' '))
print(ra, '   ("23 59 16.6")')
print(dec, '  ("-84 46 58")')
print(distance)

"""
Below are the results of your request from the Minor Planet Center's Minor Planet Ephemeris Service. Ephemerides are for the geocenter.

C/1995 O1 (Hale-Bopp)

Perturbed ephemeris below is based on elements from MPC 106342.

    CJ95O010
Date       UT      R.A. (J2000) Decl.    Delta     r     El.    Ph.   m1     Sky Motion
            h m s                                                            "/min    P.A.
2020 05 31 000000 23 59 16.6 -84 46 58  43.266  43.621  109.9   1.3  22.6    0.054   162.5
2020 06 01 000000 23 59 33.3 -84 48 12  43.265  43.625  110.1   1.3  22.6    0.054   163.4
2020 06 02 000000 23 59 49.3 -84 49 27  43.265  43.628  110.3   1.2  22.6    0.054   164.3
2020 06 03 000000 00 00 04.5 -84 50 42  43.265  43.631  110.6   1.2  22.6    0.054   165.1
2020 06 04 000000 00 00 18.9 -84 51 57  43.265  43.635  110.8   1.2  22.6    0.054   166.0
"""
