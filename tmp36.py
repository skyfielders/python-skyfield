from skyfield import api
from skyfield.api import load
from skyfield.constants import tau
from skyfield.data.mpc import COMET_URL, load_comets_dataframe

with load.open(COMET_URL) as f:
    df = load_comets_dataframe(f)

from skyfield.keplerlib import KeplerOrbit

ts = load.timescale(builtin=True)

from numpy import sqrt

df['semimajor_axis_au'] = (
    df['perihelion_distance_au'] / (1.0 - df['eccentricity'])
)

from skyfield.data.gravitational_parameters import GM_dict
from skyfield.constants import AU_KM, DAY_S

mu_km3_s2 = GM_dict[10]
mu_au3_d2 = mu_km3_s2 / (AU_KM**3.0) * (DAY_S**2.0)

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

comet = df.iloc[0:1]

k = KeplerOrbit.from_comet_dataframe(ts, comet.ix[0])

from skyfield.data.spice import inertial_frames

from math import sqrt

eph = api.load('de421.bsp')
t = ts.utc(2020, 5, 31)

k._rotation = inertial_frames['ECLIPJ2000'].T

p = eph['earth'].at(t).observe(eph['sun'] + k)
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
