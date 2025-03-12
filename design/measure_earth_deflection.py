from skyfield import api, relativity
from skyfield.api import N, W, load, wgs84

ts = load.timescale()
t = ts.utc(2025, 3, 12, 0, range(24 * 60))

eph = load('de421.bsp')
jupiter = eph['jupiter barycenter']
earth = eph['earth']
lowell = earth + wgs84.latlon(35.2029, -111.6646)

a1 = lowell.at(t).observe(jupiter).apparent()

relativity.rmasses['earth'] = 1e100  # almost no mass
a2 = lowell.at(t).observe(jupiter).apparent()

sep_mas = a2.separation_from(a1).mas()

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.set(
    title='How much Earth deflection changes Jupiter position',
    xlabel='Hours of the day 2025-03-12',
    ylabel='Angle of deflection (mas)',
)
ax.plot((t - t[0]) * 24.0, sep_mas, linestyle='--')
ax.grid()
fig.savefig('tmp.png')
