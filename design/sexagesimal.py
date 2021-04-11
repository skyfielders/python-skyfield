import math
import numpy as np
from numpy import around
from time import time

print('Which approach computes "units:minutes:sections.fraction" faster?')

places = 4
#places = 0
#places = 1

# Original approach: integers with %.

def a(value):
    t0 = time()
    power = 10 ** places
    n = int(7200.0 * power * value + 1.0) // 2
    sign = n
    n = abs(n)
    n, fraction = divmod(n, power)
    n, seconds = divmod(n, 60)

    n, minutes = divmod(n, 60)
    sign = '-' if sign < 0 else ''
    s = '%s%02dh %02dm %02d.%0*ds' % (
        sign, n, minutes, seconds, places, fraction)
    dt = time() - t0
    return dt, s

# Variant of original, just to check difference: integers with "format".

def b(value):
    t0 = time()
    power = 10 ** places
    n = m = int(7200.0 * power * value + 1.0) // 2
    n = abs(n)
    n, fraction = divmod(n, power)
    n, seconds = divmod(n, 60)

    n, minutes = divmod(n, 60)
    sign = '-' if m < 0 else ''
    s = '{0}{1:02}h {2:02}m {3:02}.{4:0{5}}s'.format(
        sign, n, minutes, seconds, fraction, places)
    dt = time() - t0
    return dt, s

# Possible new approach: floats, using round() to set the number of
# digits.  (Wow, turns out to be more expensive than I had imagined!
# The round() call is by itself the cost of either of the other
# solutions in their entirety.)

def c(value):
    t0 = time()
    value = value * 3600.0
    value = round(value, places)
    minutes, seconds = divmod(abs(value), 60.0)
    minutes = int(minutes)

    units, minutes = divmod(minutes, 60)
    sign = '-' if value < 0 else ''
    s = '{0}{1:02}h {2:02}m {3:0{4}.{5}f}s'.format(
        sign, units, minutes, seconds, places + 3, places)
    dt = time() - t0
    return dt, s

for value in 1, -1, math.pi, -0.000000001:
    value = np.array(value)
    for approach in a, b, c:
        duration, s = approach(value)
        print('{:10.8f}  {}'.format(duration, s))
    print()
