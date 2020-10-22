#!/usr/bin/env python3

import numpy as np
from skyfield.api import load

def main():
    ts = load.timescale()
    for year in range(1950, 2051, 5):
        t1 = ts.utc(year)
        t2 = ts.tt_jd(np.nextafter(t1.tt, 1e99))
        d = (t2.whole - t1.whole + t2.tt_fraction - t1.tt_fraction)
        print(year, ':', d * (24 * 60 * 60), 'seconds')

if __name__ == '__main__':
    main()
