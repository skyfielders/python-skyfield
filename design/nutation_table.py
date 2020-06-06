"""Print a table of how accurate IAU 2000B is with fewer coefficients."""

import matplotlib.pyplot as plt
from numpy import sqrt
from skyfield.api import load
from skyfield.nutationlib import iau2000a, iau2000b

ts = load.timescale()
#t = ts.utc(2000, 1, range(366 * 50))
t = ts.utc(1995, 1, range(366 * (2050 - 1995)))
print(t[0].utc_jpl())
print(t[-1].utc_jpl())

dpsi_a, depsilon_a = iau2000a(t.tt)

for i in list(range(0, 10)) + list(range(10, 71, 10)) + [77]:
    dpsi_b, depsilon_b = iau2000b(t.tt, i)
    diff1 = abs(dpsi_b - dpsi_a).max()
    diff2 = abs(depsilon_b - depsilon_a).max()
    diff_mas = sqrt(diff1 * diff1 + diff2 * diff2) * 1e-4
    if diff_mas > 1000:
        print(i, '%.3f' % (diff_mas / 1e3), 'arcseconds')
    elif diff_mas > 100:
        print(i, int(diff_mas + 0.999), 'mas')
    elif diff_mas > 10:
        print(i, '%.1f' % diff_mas, 'mas')
    else:
        print(i, '%.2f' % diff_mas, 'mas')
